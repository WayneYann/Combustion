
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
void vcomp_wdot_51_58(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
                double * restrict k_f_s, double * restrict Kc_s,
                double * restrict tc, double * restrict invT, double * restrict T);

/* Inverse molecular weights */
static const double imw[15] = {
    1.0 / 1.007970,  /*H */
    1.0 / 15.999400,  /*O */
    1.0 / 17.007370,  /*OH */
    1.0 / 2.015940,  /*H2 */
    1.0 / 31.998800,  /*O2 */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 18.015340,  /*H2O */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 30.006100,  /*NO */
    1.0 / 46.005500,  /*NO2 */
    1.0 / 44.012800,  /*N2O */
    1.0 / 15.014670,  /*NH */
    1.0 / 14.006700,  /*N */
    1.0 / 29.021370,  /*NNH */
    1.0 / 28.013400};  /*N2 */



static double fwd_A[58], fwd_beta[58], fwd_Ea[58];
static double low_A[58], low_beta[58], low_Ea[58];
static double rev_A[58], rev_beta[58], rev_Ea[58];
static double troe_a[58],troe_Ts[58], troe_Tss[58], troe_Tsss[58];
static double sri_a[58], sri_b[58], sri_c[58], sri_d[58], sri_e[58];
static double activation_units[58], prefactor_units[58], phase_units[58];
static int is_PD[58], troe_len[58], sri_len[58], nTB[58], *TBid[58];
static double *TB[58];

static double fwd_A_DEF[58], fwd_beta_DEF[58], fwd_Ea_DEF[58];
static double low_A_DEF[58], low_beta_DEF[58], low_Ea_DEF[58];
static double rev_A_DEF[58], rev_beta_DEF[58], rev_Ea_DEF[58];
static double troe_a_DEF[58],troe_Ts_DEF[58], troe_Tss_DEF[58], troe_Tsss_DEF[58];
static double sri_a_DEF[58], sri_b_DEF[58], sri_c_DEF[58], sri_d_DEF[58], sri_e_DEF[58];
static double activation_units_DEF[58], prefactor_units_DEF[58], phase_units_DEF[58];
static int is_PD_DEF[58], troe_len_DEF[58], sri_len_DEF[58], nTB_DEF[58], *TBid_DEF[58];
static double *TB_DEF[58];
static int rxn_map[58] = {9,5,10,11,12,6,0,1,7,13,14,15,8,16,17,18,19,20,21,22,23,24,25,2,26,27,28,29,30,3,31,32,33,34,4,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<58; ++i) {
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
  if (reaction_id<0 || reaction_id>=58) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=15) {
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
    for (int i=0; i<58; i++) {
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
    for (int i=0; i<58; i++) {
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
  for (int i=0; i<58; ++i) {
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
    // (0):  H + O2 <=> O + OH
    fwd_A[9]     = 3600000000000000;
    fwd_beta[9]  = -0.40999999999999998;
    fwd_Ea[9]    = 16600;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = 1e-12;
    is_PD[9] = 0;
    nTB[9] = 0;

    // (1):  H + H + M <=> H2 + M
    fwd_A[5]     = 7e+17;
    fwd_beta[5]  = -1;
    fwd_Ea[5]    = 0;
    prefactor_units[5]  = 1.0000000000000002e-12;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = 1e-12;
    is_PD[5] = 0;
    nTB[5] = 3;
    TB[5] = (double *) malloc(3 * sizeof(double));
    TBid[5] = (int *) malloc(3 * sizeof(int));
    TBid[5][0] = 14; TB[5][0] = 0; // N2
    TBid[5][1] = 6; TB[5][1] = 0; // H2O
    TBid[5][2] = 3; TB[5][2] = 0; // H2

    // (2):  H + H + N2 <=> H2 + N2
    fwd_A[10]     = 5.4e+18;
    fwd_beta[10]  = -1.3;
    fwd_Ea[10]    = 0;
    prefactor_units[10]  = 1.0000000000000002e-12;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = 1e-18;
    is_PD[10] = 0;
    nTB[10] = 0;

    // (3):  H + H + H2 <=> H2 + H2
    fwd_A[11]     = 1e+17;
    fwd_beta[11]  = -0.59999999999999998;
    fwd_Ea[11]    = 0;
    prefactor_units[11]  = 1.0000000000000002e-12;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = 1e-18;
    is_PD[11] = 0;
    nTB[11] = 0;

    // (4):  H + H + H2O <=> H2 + H2O
    fwd_A[12]     = 1e+19;
    fwd_beta[12]  = -1;
    fwd_Ea[12]    = 0;
    prefactor_units[12]  = 1.0000000000000002e-12;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = 1e-18;
    is_PD[12] = 0;
    nTB[12] = 0;

    // (5):  H + O + M <=> OH + M
    fwd_A[6]     = 62000000000000000;
    fwd_beta[6]  = -0.59999999999999998;
    fwd_Ea[6]    = 0;
    prefactor_units[6]  = 1.0000000000000002e-12;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = 1e-12;
    is_PD[6] = 0;
    nTB[6] = 1;
    TB[6] = (double *) malloc(1 * sizeof(double));
    TBid[6] = (int *) malloc(1 * sizeof(int));
    TBid[6][0] = 6; TB[6][0] = 5; // H2O

    // (6):  H + O2 (+M) <=> HO2 (+M)
    fwd_A[0]     = 1500000000000;
    fwd_beta[0]  = 0.59999999999999998;
    fwd_Ea[0]    = 0;
    low_A[0]     = 35000000000000000;
    low_beta[0]  = -0.40999999999999998;
    low_Ea[0]    = -1116;
    troe_a[0]    = 0.5;
    troe_Tsss[0] = 1.0000000000000001e-30;
    troe_Ts[0]   = 1e+30;
    troe_len[0]  = 3;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = 1e-12;
    is_PD[0] = 1;
    nTB[0] = 4;
    TB[0] = (double *) malloc(4 * sizeof(double));
    TBid[0] = (int *) malloc(4 * sizeof(int));
    TBid[0][0] = 14; TB[0][0] = 0; // N2
    TBid[0][1] = 6; TB[0][1] = 11; // H2O
    TBid[0][2] = 3; TB[0][2] = 2; // H2
    TBid[0][3] = 4; TB[0][3] = 0.78000000000000003; // O2

    // (7):  H + O2 (+N2) <=> HO2 (+N2)
    fwd_A[1]     = 1500000000000;
    fwd_beta[1]  = 0.59999999999999998;
    fwd_Ea[1]    = 0;
    low_A[1]     = 6.37e+20;
    low_beta[1]  = -1.72;
    low_Ea[1]    = 520;
    troe_a[1]    = 0.80000000000000004;
    troe_Tsss[1] = 1.0000000000000001e-30;
    troe_Ts[1]   = 1e+30;
    troe_len[1]  = 3;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = 1e-12;
    is_PD[1] = 1;
    nTB[1] = 0;
    TB[1] = (double *) malloc(0 * sizeof(double));
    TBid[1] = (int *) malloc(0 * sizeof(int));

    // (8):  O + O + M <=> O2 + M
    fwd_A[7]     = 19000000000000;
    fwd_beta[7]  = 0;
    fwd_Ea[7]    = -1788;
    prefactor_units[7]  = 1.0000000000000002e-12;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = 1e-12;
    is_PD[7] = 0;
    nTB[7] = 3;
    TB[7] = (double *) malloc(3 * sizeof(double));
    TBid[7] = (int *) malloc(3 * sizeof(int));
    TBid[7][0] = 14; TB[7][0] = 1.5; // N2
    TBid[7][1] = 4; TB[7][1] = 1.5; // O2
    TBid[7][2] = 6; TB[7][2] = 10; // H2O

    // (9):  O + H2 <=> OH + H
    fwd_A[13]     = 3800000000000;
    fwd_beta[13]  = 0;
    fwd_Ea[13]    = 7948;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = 1e-12;
    is_PD[13] = 0;
    nTB[13] = 0;

    // (10):  O + H2 <=> OH + H
    fwd_A[14]     = 880000000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = 19175;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = 1e-12;
    is_PD[14] = 0;
    nTB[14] = 0;

    // (11):  OH + OH <=> O + H2O
    fwd_A[15]     = 4300;
    fwd_beta[15]  = 2.7000000000000002;
    fwd_Ea[15]    = -1822;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = 1e-12;
    is_PD[15] = 0;
    nTB[15] = 0;

    // (12):  OH + H + M <=> H2O + M
    fwd_A[8]     = 4.4999999999999998e+22;
    fwd_beta[8]  = -2;
    fwd_Ea[8]    = 0;
    prefactor_units[8]  = 1.0000000000000002e-12;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = 1e-12;
    is_PD[8] = 0;
    nTB[8] = 2;
    TB[8] = (double *) malloc(2 * sizeof(double));
    TBid[8] = (int *) malloc(2 * sizeof(int));
    TBid[8][0] = 3; TB[8][0] = 0.72999999999999998; // H2
    TBid[8][1] = 6; TB[8][1] = 12; // H2O

    // (13):  OH + H2 <=> H + H2O
    fwd_A[16]     = 210000000;
    fwd_beta[16]  = 1.52;
    fwd_Ea[16]    = 3449;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = 1e-12;
    is_PD[16] = 0;
    nTB[16] = 0;

    // (14):  H2 + O2 <=> HO2 + H
    fwd_A[17]     = 740000;
    fwd_beta[17]  = 2.4329999999999998;
    fwd_Ea[17]    = 53502;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = 1e-12;
    is_PD[17] = 0;
    nTB[17] = 0;

    // (15):  HO2 + H <=> OH + OH
    fwd_A[18]     = 84000000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 400;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = 1e-12;
    is_PD[18] = 0;
    nTB[18] = 0;

    // (16):  HO2 + H <=> H2O + O
    fwd_A[19]     = 1400000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = 0;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = 1e-12;
    is_PD[19] = 0;
    nTB[19] = 0;

    // (17):  HO2 + O <=> OH + O2
    fwd_A[20]     = 16000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = -445;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = 1e-12;
    is_PD[20] = 0;
    nTB[20] = 0;

    // (18):  HO2 + OH <=> H2O + O2
    fwd_A[21]     = 3.6e+21;
    fwd_beta[21]  = -2.1000000000000001;
    fwd_Ea[21]    = 9000;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = 1e-12;
    is_PD[21] = 0;
    nTB[21] = 0;

    // (19):  HO2 + OH <=> H2O + O2
    fwd_A[22]     = 2000000000000000;
    fwd_beta[22]  = -0.59999999999999998;
    fwd_Ea[22]    = 0;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = 1e-12;
    is_PD[22] = 0;
    nTB[22] = 0;

    // (20):  HO2 + OH <=> H2O + O2
    fwd_A[23]     = -2.1999999999999999e+96;
    fwd_beta[23]  = -24;
    fwd_Ea[23]    = 49000;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = 1e-12;
    is_PD[23] = 0;
    nTB[23] = 0;

    // (21):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[24]     = 190000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = -1408;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = 1e-12;
    is_PD[24] = 0;
    nTB[24] = 0;

    // (22):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[25]     = 100000000000000;
    fwd_beta[25]  = 0;
    fwd_Ea[25]    = 11034;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = 1e-12;
    is_PD[25] = 0;
    nTB[25] = 0;

    // (23):  H2O2 (+M) <=> OH + OH (+M)
    fwd_A[2]     = 400000000000;
    fwd_beta[2]  = 0;
    fwd_Ea[2]    = 37137;
    low_A[2]     = 22910000000000000;
    low_beta[2]  = 0;
    low_Ea[2]    = 43638;
    troe_a[2]    = 0.5;
    troe_Tsss[2] = 1.0000000000000001e-30;
    troe_Ts[2]   = 1e+30;
    troe_Tss[2]  = 1e+30;
    troe_len[2]  = 4;
    prefactor_units[2]  = 1;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = 1e-6;
    is_PD[2] = 1;
    nTB[2] = 2;
    TB[2] = (double *) malloc(2 * sizeof(double));
    TBid[2] = (int *) malloc(2 * sizeof(int));
    TBid[2][0] = 6; TB[2][0] = 12; // H2O
    TBid[2][1] = 3; TB[2][1] = 2.5; // H2

    // (24):  H2O2 + H <=> H2O + OH
    fwd_A[26]     = 10000000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 3580;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = 1e-12;
    is_PD[26] = 0;
    nTB[26] = 0;

    // (25):  H2O2 + H <=> HO2 + H2
    fwd_A[27]     = 1700000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = 3760;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = 1e-12;
    is_PD[27] = 0;
    nTB[27] = 0;

    // (26):  H2O2 + O <=> HO2 + OH
    fwd_A[28]     = 9600000;
    fwd_beta[28]  = 2;
    fwd_Ea[28]    = 3970;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = 1e-12;
    is_PD[28] = 0;
    nTB[28] = 0;

    // (27):  H2O2 + OH <=> H2O + HO2
    fwd_A[29]     = 1900000000000;
    fwd_beta[29]  = 0;
    fwd_Ea[29]    = 427;
    prefactor_units[29]  = 1.0000000000000002e-06;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = 1e-12;
    is_PD[29] = 0;
    nTB[29] = 0;

    // (28):  H2O2 + OH <=> H2O + HO2
    fwd_A[30]     = 1.6e+18;
    fwd_beta[30]  = 0;
    fwd_Ea[30]    = 29410;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = 1e-12;
    is_PD[30] = 0;
    nTB[30] = 0;

    // (29):  NO + O (+M) <=> NO2 (+M)
    fwd_A[3]     = 1300000000000000;
    fwd_beta[3]  = -0.75;
    fwd_Ea[3]    = 0;
    low_A[3]     = 4.7199999999999999e+24;
    low_beta[3]  = -2.8700000000000001;
    low_Ea[3]    = 1550;
    troe_a[3]    = 0.88;
    troe_Tsss[3] = 1000;
    troe_Ts[3]   = 10000;
    troe_Tss[3]  = 1e+30;
    troe_len[3]  = 4;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = 1e-12;
    is_PD[3] = 1;
    nTB[3] = 0;
    TB[3] = (double *) malloc(0 * sizeof(double));
    TBid[3] = (int *) malloc(0 * sizeof(int));

    // (30):  NO + HO2 <=> NO2 + OH
    fwd_A[31]     = 2100000000000;
    fwd_beta[31]  = 0;
    fwd_Ea[31]    = -497;
    prefactor_units[31]  = 1.0000000000000002e-06;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = 1e-12;
    is_PD[31] = 0;
    nTB[31] = 0;

    // (31):  NO2 + H <=> NO + OH
    fwd_A[32]     = 130000000000000;
    fwd_beta[32]  = 0;
    fwd_Ea[32]    = 362;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = 1e-12;
    is_PD[32] = 0;
    nTB[32] = 0;

    // (32):  NO2 + O <=> NO + O2
    fwd_A[33]     = 110000000000000;
    fwd_beta[33]  = -0.52000000000000002;
    fwd_Ea[33]    = 0;
    prefactor_units[33]  = 1.0000000000000002e-06;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = 1e-12;
    is_PD[33] = 0;
    nTB[33] = 0;

    // (33):  NO2 + NO2 <=> NO + NO + O2
    fwd_A[34]     = 4500000000000;
    fwd_beta[34]  = 0;
    fwd_Ea[34]    = 27599;
    prefactor_units[34]  = 1.0000000000000002e-06;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = 1e-12;
    is_PD[34] = 0;
    nTB[34] = 0;

    // (34):  N2O (+M) <=> N2 + O (+M)
    fwd_A[4]     = 1300000000000;
    fwd_beta[4]  = 0;
    fwd_Ea[4]    = 62570;
    low_A[4]     = 400000000000000;
    low_beta[4]  = 0;
    low_Ea[4]    = 56600;
    prefactor_units[4]  = 1;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = 1e-6;
    is_PD[4] = 1;
    nTB[4] = 3;
    TB[4] = (double *) malloc(3 * sizeof(double));
    TBid[4] = (int *) malloc(3 * sizeof(int));
    TBid[4][0] = 14; TB[4][0] = 1.7; // N2
    TBid[4][1] = 4; TB[4][1] = 1.3999999999999999; // O2
    TBid[4][2] = 6; TB[4][2] = 12; // H2O

    // (35):  N2O + H <=> N2 + OH
    fwd_A[35]     = 33000000000;
    fwd_beta[35]  = 0;
    fwd_Ea[35]    = 4729;
    prefactor_units[35]  = 1.0000000000000002e-06;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = 1e-12;
    is_PD[35] = 0;
    nTB[35] = 0;

    // (36):  N2O + H <=> N2 + OH
    fwd_A[36]     = 440000000000000;
    fwd_beta[36]  = 0;
    fwd_Ea[36]    = 19254;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = 1e-12;
    is_PD[36] = 0;
    nTB[36] = 0;

    // (37):  N2O + O <=> NO + NO
    fwd_A[37]     = 92000000000000;
    fwd_beta[37]  = 0;
    fwd_Ea[37]    = 27679;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = 1e-12;
    is_PD[37] = 0;
    nTB[37] = 0;

    // (38):  N2O + O <=> N2 + O2
    fwd_A[38]     = 3700000000000;
    fwd_beta[38]  = 0;
    fwd_Ea[38]    = 15936;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = 1e-12;
    is_PD[38] = 0;
    nTB[38] = 0;

    // (39):  NH + H <=> N + H2
    fwd_A[39]     = 30000000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 0;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = 1e-12;
    is_PD[39] = 0;
    nTB[39] = 0;

    // (40):  NH + O <=> NO + H
    fwd_A[40]     = 92000000000000;
    fwd_beta[40]  = 0;
    fwd_Ea[40]    = 0;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = 1e-12;
    is_PD[40] = 0;
    nTB[40] = 0;

    // (41):  NH + OH <=> N + H2O
    fwd_A[41]     = 500000000000;
    fwd_beta[41]  = 0.5;
    fwd_Ea[41]    = 2000;
    prefactor_units[41]  = 1.0000000000000002e-06;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = 1e-12;
    is_PD[41] = 0;
    nTB[41] = 0;

    // (42):  NH + O2 <=> NO + OH
    fwd_A[42]     = 1300000;
    fwd_beta[42]  = 1.5;
    fwd_Ea[42]    = 100;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = 1e-12;
    is_PD[42] = 0;
    nTB[42] = 0;

    // (43):  NH + NO <=> N2O + H
    fwd_A[43]     = 290000000000000;
    fwd_beta[43]  = -0.40000000000000002;
    fwd_Ea[43]    = 0;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = 1e-12;
    is_PD[43] = 0;
    nTB[43] = 0;

    // (44):  NH + NO <=> N2O + H
    fwd_A[44]     = -22000000000000;
    fwd_beta[44]  = -0.23000000000000001;
    fwd_Ea[44]    = 0;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = 1e-12;
    is_PD[44] = 0;
    nTB[44] = 0;

    // (45):  NH + NO <=> N2 + OH
    fwd_A[45]     = 22000000000000;
    fwd_beta[45]  = -0.23000000000000001;
    fwd_Ea[45]    = 0;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = 1e-12;
    is_PD[45] = 0;
    nTB[45] = 0;

    // (46):  NH + NO2 <=> N2O + OH
    fwd_A[46]     = 10000000000000;
    fwd_beta[46]  = 0;
    fwd_Ea[46]    = 0;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = 1e-12;
    is_PD[46] = 0;
    nTB[46] = 0;

    // (47):  N + OH <=> NO + H
    fwd_A[47]     = 38000000000000;
    fwd_beta[47]  = 0;
    fwd_Ea[47]    = 0;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = 1e-12;
    is_PD[47] = 0;
    nTB[47] = 0;

    // (48):  N + O2 <=> NO + O
    fwd_A[48]     = 6400000000;
    fwd_beta[48]  = 1;
    fwd_Ea[48]    = 6280;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = 1e-12;
    is_PD[48] = 0;
    nTB[48] = 0;

    // (49):  N + NO <=> N2 + O
    fwd_A[49]     = 21000000000000;
    fwd_beta[49]  = 0;
    fwd_Ea[49]    = 0;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = 1e-12;
    is_PD[49] = 0;
    nTB[49] = 0;

    // (50):  NNH <=> N2 + H
    fwd_A[50]     = 65000000;
    fwd_beta[50]  = 0;
    fwd_Ea[50]    = 0;
    prefactor_units[50]  = 1;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = 1e-6;
    is_PD[50] = 0;
    nTB[50] = 0;

    // (51):  NNH + H <=> N2 + H2
    fwd_A[51]     = 100000000000000;
    fwd_beta[51]  = 0;
    fwd_Ea[51]    = 0;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = 1e-12;
    is_PD[51] = 0;
    nTB[51] = 0;

    // (52):  NNH + O <=> N2O + H
    fwd_A[52]     = 100000000000000;
    fwd_beta[52]  = 0;
    fwd_Ea[52]    = 0;
    prefactor_units[52]  = 1.0000000000000002e-06;
    activation_units[52] = 0.50321666580471969;
    phase_units[52]      = 1e-12;
    is_PD[52] = 0;
    nTB[52] = 0;

    // (53):  NNH + O <=> N2 + OH
    fwd_A[53]     = 80000000000000;
    fwd_beta[53]  = 0;
    fwd_Ea[53]    = 0;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = 1e-12;
    is_PD[53] = 0;
    nTB[53] = 0;

    // (54):  NNH + O <=> NH + NO
    fwd_A[54]     = 50000000000000;
    fwd_beta[54]  = 0;
    fwd_Ea[54]    = 0;
    prefactor_units[54]  = 1.0000000000000002e-06;
    activation_units[54] = 0.50321666580471969;
    phase_units[54]      = 1e-12;
    is_PD[54] = 0;
    nTB[54] = 0;

    // (55):  NNH + OH <=> N2 + H2O
    fwd_A[55]     = 50000000000000;
    fwd_beta[55]  = 0;
    fwd_Ea[55]    = 0;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = 1e-12;
    is_PD[55] = 0;
    nTB[55] = 0;

    // (56):  NNH + O2 <=> N2 + HO2
    fwd_A[56]     = 200000000000000;
    fwd_beta[56]  = 0;
    fwd_Ea[56]    = 0;
    prefactor_units[56]  = 1.0000000000000002e-06;
    activation_units[56] = 0.50321666580471969;
    phase_units[56]      = 1e-12;
    is_PD[56] = 0;
    nTB[56] = 0;

    // (57):  NNH + O2 <=> N2 + H + O2
    fwd_A[57]     = 50000000000000;
    fwd_beta[57]  = 0;
    fwd_Ea[57]    = 0;
    prefactor_units[57]  = 1.0000000000000002e-06;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = 1e-12;
    is_PD[57] = 0;
    nTB[57] = 0;

    SetAllDefaults();
}



/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 3;
    *kk = 15;
    *ii = 58;
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
    for (i=0; i<lenkname*15; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = 'H';
    kname[ 2*lenkname + 2 ] = ' ';

    /* H2  */
    kname[ 3*lenkname + 0 ] = 'H';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* O2  */
    kname[ 4*lenkname + 0 ] = 'O';
    kname[ 4*lenkname + 1 ] = '2';
    kname[ 4*lenkname + 2 ] = ' ';

    /* HO2  */
    kname[ 5*lenkname + 0 ] = 'H';
    kname[ 5*lenkname + 1 ] = 'O';
    kname[ 5*lenkname + 2 ] = '2';
    kname[ 5*lenkname + 3 ] = ' ';

    /* H2O  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = '2';
    kname[ 6*lenkname + 2 ] = 'O';
    kname[ 6*lenkname + 3 ] = ' ';

    /* H2O2  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = '2';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = '2';
    kname[ 7*lenkname + 4 ] = ' ';

    /* NO  */
    kname[ 8*lenkname + 0 ] = 'N';
    kname[ 8*lenkname + 1 ] = 'O';
    kname[ 8*lenkname + 2 ] = ' ';

    /* NO2  */
    kname[ 9*lenkname + 0 ] = 'N';
    kname[ 9*lenkname + 1 ] = 'O';
    kname[ 9*lenkname + 2 ] = '2';
    kname[ 9*lenkname + 3 ] = ' ';

    /* N2O  */
    kname[ 10*lenkname + 0 ] = 'N';
    kname[ 10*lenkname + 1 ] = '2';
    kname[ 10*lenkname + 2 ] = 'O';
    kname[ 10*lenkname + 3 ] = ' ';

    /* NH  */
    kname[ 11*lenkname + 0 ] = 'N';
    kname[ 11*lenkname + 1 ] = 'H';
    kname[ 11*lenkname + 2 ] = ' ';

    /* N  */
    kname[ 12*lenkname + 0 ] = 'N';
    kname[ 12*lenkname + 1 ] = ' ';

    /* NNH  */
    kname[ 13*lenkname + 0 ] = 'N';
    kname[ 13*lenkname + 1 ] = 'N';
    kname[ 13*lenkname + 2 ] = 'H';
    kname[ 13*lenkname + 3 ] = ' ';

    /* N2  */
    kname[ 14*lenkname + 0 ] = 'N';
    kname[ 14*lenkname + 1 ] = '2';
    kname[ 14*lenkname + 2 ] = ' ';

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
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*O */
    YOW += y[2]*imw[2]; /*OH */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*O2 */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*NO */
    YOW += y[9]*imw[9]; /*NO2 */
    YOW += y[10]*imw[10]; /*N2O */
    YOW += y[11]*imw[11]; /*NH */
    YOW += y[12]*imw[12]; /*N */
    YOW += y[13]*imw[13]; /*NNH */
    YOW += y[14]*imw[14]; /*N2 */
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

    for (int n=0; n<15; n++) {
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
    W += c[0]*1.007970; /*H */
    W += c[1]*15.999400; /*O */
    W += c[2]*17.007370; /*OH */
    W += c[3]*2.015940; /*H2 */
    W += c[4]*31.998800; /*O2 */
    W += c[5]*33.006770; /*HO2 */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*30.006100; /*NO */
    W += c[9]*46.005500; /*NO2 */
    W += c[10]*44.012800; /*N2O */
    W += c[11]*15.014670; /*NH */
    W += c[12]*14.006700; /*N */
    W += c[13]*29.021370; /*NNH */
    W += c[14]*28.013400; /*N2 */

    for (id = 0; id < 15; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double YOW = 0;
    double tmp[15];

    for (int i = 0; i < 15; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 15; i++)
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
    W += c[0]*1.007970; /*H */
    W += c[1]*15.999400; /*O */
    W += c[2]*17.007370; /*OH */
    W += c[3]*2.015940; /*H2 */
    W += c[4]*31.998800; /*O2 */
    W += c[5]*33.006770; /*HO2 */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*30.006100; /*NO */
    W += c[9]*46.005500; /*NO2 */
    W += c[10]*44.012800; /*N2O */
    W += c[11]*15.014670; /*NH */
    W += c[12]*14.006700; /*N */
    W += c[13]*29.021370; /*NNH */
    W += c[14]*28.013400; /*N2 */

    for (id = 0; id < 15; ++id) {
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
    double tmp[15];

    for (int i = 0; i < 15; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 15; i++)
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
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
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
    W += c[0]*1.007970; /*H */
    W += c[1]*15.999400; /*O */
    W += c[2]*17.007370; /*OH */
    W += c[3]*2.015940; /*H2 */
    W += c[4]*31.998800; /*O2 */
    W += c[5]*33.006770; /*HO2 */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*30.006100; /*NO */
    W += c[9]*46.005500; /*NO2 */
    W += c[10]*44.012800; /*N2O */
    W += c[11]*15.014670; /*NH */
    W += c[12]*14.006700; /*N */
    W += c[13]*29.021370; /*NNH */
    W += c[14]*28.013400; /*N2 */

    for (id = 0; id < 15; ++id) {
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
    double tmp[15];

    for (int i = 0; i < 15; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 15; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 15; i++)
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

    for (int n=0; n<15; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<15; n++) {
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
    for (int i = 0; i < 15; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 15; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 15; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int i = 0; i < 15; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*1.007970*XWinv; 
    y[1] = x[1]*15.999400*XWinv; 
    y[2] = x[2]*17.007370*XWinv; 
    y[3] = x[3]*2.015940*XWinv; 
    y[4] = x[4]*31.998800*XWinv; 
    y[5] = x[5]*33.006770*XWinv; 
    y[6] = x[6]*18.015340*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*30.006100*XWinv; 
    y[9] = x[9]*46.005500*XWinv; 
    y[10] = x[10]*44.012800*XWinv; 
    y[11] = x[11]*15.014670*XWinv; 
    y[12] = x[12]*14.006700*XWinv; 
    y[13] = x[13]*29.021370*XWinv; 
    y[14] = x[14]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 15; ++id) {
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
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 15; ++id) {
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
    for (id = 0; id < 15; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 15; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*1.007970; /*H */
    CW += c[1]*15.999400; /*O */
    CW += c[2]*17.007370; /*OH */
    CW += c[3]*2.015940; /*H2 */
    CW += c[4]*31.998800; /*O2 */
    CW += c[5]*33.006770; /*HO2 */
    CW += c[6]*18.015340; /*H2O */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*30.006100; /*NO */
    CW += c[9]*46.005500; /*NO2 */
    CW += c[10]*44.012800; /*N2O */
    CW += c[11]*15.014670; /*NH */
    CW += c[12]*14.006700; /*N */
    CW += c[13]*29.021370; /*NNH */
    CW += c[14]*28.013400; /*N2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*1.007970*CWinv; 
    y[1] = c[1]*15.999400*CWinv; 
    y[2] = c[2]*17.007370*CWinv; 
    y[3] = c[3]*2.015940*CWinv; 
    y[4] = c[4]*31.998800*CWinv; 
    y[5] = c[5]*33.006770*CWinv; 
    y[6] = c[6]*18.015340*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*30.006100*CWinv; 
    y[9] = c[9]*46.005500*CWinv; 
    y[10] = c[10]*44.012800*CWinv; 
    y[11] = c[11]*15.014670*CWinv; 
    y[12] = c[12]*14.006700*CWinv; 
    y[13] = c[13]*29.021370*CWinv; 
    y[14] = c[14]*28.013400*CWinv; 

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
    for (id = 0; id < 15; ++id) {
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
    for (id = 0; id < 15; ++id) {
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
    for (id = 0; id < 15; ++id) {
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
    for (id = 0; id < 15; ++id) {
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
    for (id = 0; id < 15; ++id) {
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
    for (id = 0; id < 15; ++id) {
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
    for (id = 0; id < 15; ++id) {
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
    cvms[0] *= 8.248767324424338e+07; /*H */
    cvms[1] *= 5.196763628636074e+06; /*O */
    cvms[2] *= 4.888768810227566e+06; /*OH */
    cvms[3] *= 4.124383662212169e+07; /*H2 */
    cvms[4] *= 2.598381814318037e+06; /*O2 */
    cvms[5] *= 2.519031701678171e+06; /*HO2 */
    cvms[6] *= 4.615239012974499e+06; /*H2O */
    cvms[7] *= 2.444384405113783e+06; /*H2O2 */
    cvms[8] *= 2.770939908885194e+06; /*NO */
    cvms[9] *= 1.807286085359359e+06; /*NO2 */
    cvms[10] *= 1.889111803838883e+06; /*N2O */
    cvms[11] *= 5.537590902763763e+06; /*NH */
    cvms[12] *= 5.936094868884177e+06; /*N */
    cvms[13] *= 2.864961233739138e+06; /*NNH */
    cvms[14] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 8.248767324424338e+07; /*H */
    cpms[1] *= 5.196763628636074e+06; /*O */
    cpms[2] *= 4.888768810227566e+06; /*OH */
    cpms[3] *= 4.124383662212169e+07; /*H2 */
    cpms[4] *= 2.598381814318037e+06; /*O2 */
    cpms[5] *= 2.519031701678171e+06; /*HO2 */
    cpms[6] *= 4.615239012974499e+06; /*H2O */
    cpms[7] *= 2.444384405113783e+06; /*H2O2 */
    cpms[8] *= 2.770939908885194e+06; /*NO */
    cpms[9] *= 1.807286085359359e+06; /*NO2 */
    cpms[10] *= 1.889111803838883e+06; /*N2O */
    cpms[11] *= 5.537590902763763e+06; /*NH */
    cpms[12] *= 5.936094868884177e+06; /*N */
    cpms[13] *= 2.864961233739138e+06; /*NNH */
    cpms[14] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 15; i++)
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
    for (int i = 0; i < 15; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tc[5], h[15];

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
    }

    for (int n=0; n<15; n++) {
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
    for (int i = 0; i < 15; i++)
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
    for (int i = 0; i < 15; i++)
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
    sms[0] *= 8.248767324424338e+07; /*H */
    sms[1] *= 5.196763628636074e+06; /*O */
    sms[2] *= 4.888768810227566e+06; /*OH */
    sms[3] *= 4.124383662212169e+07; /*H2 */
    sms[4] *= 2.598381814318037e+06; /*O2 */
    sms[5] *= 2.519031701678171e+06; /*HO2 */
    sms[6] *= 4.615239012974499e+06; /*H2O */
    sms[7] *= 2.444384405113783e+06; /*H2O2 */
    sms[8] *= 2.770939908885194e+06; /*NO */
    sms[9] *= 1.807286085359359e+06; /*NO2 */
    sms[10] *= 1.889111803838883e+06; /*N2O */
    sms[11] *= 5.537590902763763e+06; /*NH */
    sms[12] *= 5.936094868884177e+06; /*N */
    sms[13] *= 2.864961233739138e+06; /*NNH */
    sms[14] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[15]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 15; ++id) {
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
    double cpor[15], tresult[15]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 15; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 15; i++)
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
    double cvor[15]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 15; ++id) {
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
    double cvor[15]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H */
    result += cvor[1]*y[1]*imw[1]; /*O */
    result += cvor[2]*y[2]*imw[2]; /*OH */
    result += cvor[3]*y[3]*imw[3]; /*H2 */
    result += cvor[4]*y[4]*imw[4]; /*O2 */
    result += cvor[5]*y[5]*imw[5]; /*HO2 */
    result += cvor[6]*y[6]*imw[6]; /*H2O */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*NO */
    result += cvor[9]*y[9]*imw[9]; /*NO2 */
    result += cvor[10]*y[10]*imw[10]; /*N2O */
    result += cvor[11]*y[11]*imw[11]; /*NH */
    result += cvor[12]*y[12]*imw[12]; /*N */
    result += cvor[13]*y[13]*imw[13]; /*NNH */
    result += cvor[14]*y[14]*imw[14]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[15]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 15; ++id) {
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
    double hml[15], tmp[15]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 15; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 15; ++id) {
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
    double uml[15]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 15; ++id) {
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
    double ums[15]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H */
    result += y[1]*ums[1]*imw[1]; /*O */
    result += y[2]*ums[2]*imw[2]; /*OH */
    result += y[3]*ums[3]*imw[3]; /*H2 */
    result += y[4]*ums[4]*imw[4]; /*O2 */
    result += y[5]*ums[5]*imw[5]; /*HO2 */
    result += y[6]*ums[6]*imw[6]; /*H2O */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*NO */
    result += y[9]*ums[9]*imw[9]; /*NO2 */
    result += y[10]*ums[10]*imw[10]; /*N2O */
    result += y[11]*ums[11]*imw[11]; /*NH */
    result += y[12]*ums[12]*imw[12]; /*N */
    result += y[13]*ums[13]*imw[13]; /*NNH */
    result += y[14]*ums[14]*imw[14]; /*N2 */

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
    double sor[15]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 15; ++id) {
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
    double sor[15]; /* temporary storage */
    double x[15]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*O */
    YOW += y[2]*imw[2]; /*OH */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*O2 */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*NO */
    YOW += y[9]*imw[9]; /*NO2 */
    YOW += y[10]*imw[10]; /*N2O */
    YOW += y[11]*imw[11]; /*NH */
    YOW += y[12]*imw[12]; /*N */
    YOW += y[13]*imw[13]; /*NNH */
    YOW += y[14]*imw[14]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(15.999400*YOW); 
    x[2] = y[2]/(17.007370*YOW); 
    x[3] = y[3]/(2.015940*YOW); 
    x[4] = y[4]/(31.998800*YOW); 
    x[5] = y[5]/(33.006770*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(30.006100*YOW); 
    x[9] = y[9]/(46.005500*YOW); 
    x[10] = y[10]/(44.012800*YOW); 
    x[11] = y[11]/(15.014670*YOW); 
    x[12] = y[12]/(14.006700*YOW); 
    x[13] = y[13]/(29.021370*YOW); 
    x[14] = y[14]/(28.013400*YOW); 
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
    double gort[15]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 15; ++id) {
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
    double gort[15]; /* temporary storage */
    double x[15]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*O */
    YOW += y[2]*imw[2]; /*OH */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*O2 */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*NO */
    YOW += y[9]*imw[9]; /*NO2 */
    YOW += y[10]*imw[10]; /*N2O */
    YOW += y[11]*imw[11]; /*NH */
    YOW += y[12]*imw[12]; /*N */
    YOW += y[13]*imw[13]; /*NNH */
    YOW += y[14]*imw[14]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(15.999400*YOW); 
    x[2] = y[2]/(17.007370*YOW); 
    x[3] = y[3]/(2.015940*YOW); 
    x[4] = y[4]/(31.998800*YOW); 
    x[5] = y[5]/(33.006770*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(30.006100*YOW); 
    x[9] = y[9]/(46.005500*YOW); 
    x[10] = y[10]/(44.012800*YOW); 
    x[11] = y[11]/(15.014670*YOW); 
    x[12] = y[12]/(14.006700*YOW); 
    x[13] = y[13]/(29.021370*YOW); 
    x[14] = y[14]/(28.013400*YOW); 
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
    double aort[15]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 15; ++id) {
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
    double aort[15]; /* temporary storage */
    double x[15]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*O */
    YOW += y[2]*imw[2]; /*OH */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*O2 */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*NO */
    YOW += y[9]*imw[9]; /*NO2 */
    YOW += y[10]*imw[10]; /*N2O */
    YOW += y[11]*imw[11]; /*NH */
    YOW += y[12]*imw[12]; /*N */
    YOW += y[13]*imw[13]; /*NNH */
    YOW += y[14]*imw[14]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(15.999400*YOW); 
    x[2] = y[2]/(17.007370*YOW); 
    x[3] = y[3]/(2.015940*YOW); 
    x[4] = y[4]/(31.998800*YOW); 
    x[5] = y[5]/(33.006770*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(30.006100*YOW); 
    x[9] = y[9]/(46.005500*YOW); 
    x[10] = y[10]/(44.012800*YOW); 
    x[11] = y[11]/(15.014670*YOW); 
    x[12] = y[12]/(14.006700*YOW); 
    x[13] = y[13]/(29.021370*YOW); 
    x[14] = y[14]/(28.013400*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 15; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*O */
    YOW += y[2]*imw[2]; /*OH */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*O2 */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*NO */
    YOW += y[9]*imw[9]; /*NO2 */
    YOW += y[10]*imw[10]; /*N2O */
    YOW += y[11]*imw[11]; /*NH */
    YOW += y[12]*imw[12]; /*N */
    YOW += y[13]*imw[13]; /*NNH */
    YOW += y[14]*imw[14]; /*N2 */
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

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
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

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[15*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<15; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<15*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 15; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 58; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 58; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*O */
    YOW += y[2]*imw[2]; /*OH */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*O2 */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*NO */
    YOW += y[9]*imw[9]; /*NO2 */
    YOW += y[10]*imw[10]; /*N2O */
    YOW += y[11]*imw[11]; /*NH */
    YOW += y[12]*imw[12]; /*N */
    YOW += y[13]*imw[13]; /*NNH */
    YOW += y[14]*imw[14]; /*N2 */
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

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 58; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 58; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
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

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 58; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 58; ++id) {
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
    for (id = 0; id < 15 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 0 * kd + 0 ] += -1 ;
    nuki[ 4 * kd + 0 ] += -1 ;
    nuki[ 5 * kd + 0 ] += +1 ;

    /*reaction 2: H + O2 (+N2) <=> HO2 (+N2) */
    nuki[ 0 * kd + 1 ] += -1 ;
    nuki[ 4 * kd + 1 ] += -1 ;
    nuki[ 5 * kd + 1 ] += +1 ;

    /*reaction 3: H2O2 (+M) <=> OH + OH (+M) */
    nuki[ 7 * kd + 2 ] += -1 ;
    nuki[ 2 * kd + 2 ] += +1 ;
    nuki[ 2 * kd + 2 ] += +1 ;

    /*reaction 4: NO + O (+M) <=> NO2 (+M) */
    nuki[ 8 * kd + 3 ] += -1 ;
    nuki[ 1 * kd + 3 ] += -1 ;
    nuki[ 9 * kd + 3 ] += +1 ;

    /*reaction 5: N2O (+M) <=> N2 + O (+M) */
    nuki[ 10 * kd + 4 ] += -1 ;
    nuki[ 14 * kd + 4 ] += +1 ;
    nuki[ 1 * kd + 4 ] += +1 ;

    /*reaction 6: H + H + M <=> H2 + M */
    nuki[ 0 * kd + 5 ] += -1 ;
    nuki[ 0 * kd + 5 ] += -1 ;
    nuki[ 3 * kd + 5 ] += +1 ;

    /*reaction 7: H + O + M <=> OH + M */
    nuki[ 0 * kd + 6 ] += -1 ;
    nuki[ 1 * kd + 6 ] += -1 ;
    nuki[ 2 * kd + 6 ] += +1 ;

    /*reaction 8: O + O + M <=> O2 + M */
    nuki[ 1 * kd + 7 ] += -1 ;
    nuki[ 1 * kd + 7 ] += -1 ;
    nuki[ 4 * kd + 7 ] += +1 ;

    /*reaction 9: OH + H + M <=> H2O + M */
    nuki[ 2 * kd + 8 ] += -1 ;
    nuki[ 0 * kd + 8 ] += -1 ;
    nuki[ 6 * kd + 8 ] += +1 ;

    /*reaction 10: H + O2 <=> O + OH */
    nuki[ 0 * kd + 9 ] += -1 ;
    nuki[ 4 * kd + 9 ] += -1 ;
    nuki[ 1 * kd + 9 ] += +1 ;
    nuki[ 2 * kd + 9 ] += +1 ;

    /*reaction 11: H + H + N2 <=> H2 + N2 */
    nuki[ 0 * kd + 10 ] += -1 ;
    nuki[ 0 * kd + 10 ] += -1 ;
    nuki[ 14 * kd + 10 ] += -1 ;
    nuki[ 3 * kd + 10 ] += +1 ;
    nuki[ 14 * kd + 10 ] += +1 ;

    /*reaction 12: H + H + H2 <=> H2 + H2 */
    nuki[ 0 * kd + 11 ] += -1 ;
    nuki[ 0 * kd + 11 ] += -1 ;
    nuki[ 3 * kd + 11 ] += -1 ;
    nuki[ 3 * kd + 11 ] += +1 ;
    nuki[ 3 * kd + 11 ] += +1 ;

    /*reaction 13: H + H + H2O <=> H2 + H2O */
    nuki[ 0 * kd + 12 ] += -1 ;
    nuki[ 0 * kd + 12 ] += -1 ;
    nuki[ 6 * kd + 12 ] += -1 ;
    nuki[ 3 * kd + 12 ] += +1 ;
    nuki[ 6 * kd + 12 ] += +1 ;

    /*reaction 14: O + H2 <=> OH + H */
    nuki[ 1 * kd + 13 ] += -1 ;
    nuki[ 3 * kd + 13 ] += -1 ;
    nuki[ 2 * kd + 13 ] += +1 ;
    nuki[ 0 * kd + 13 ] += +1 ;

    /*reaction 15: O + H2 <=> OH + H */
    nuki[ 1 * kd + 14 ] += -1 ;
    nuki[ 3 * kd + 14 ] += -1 ;
    nuki[ 2 * kd + 14 ] += +1 ;
    nuki[ 0 * kd + 14 ] += +1 ;

    /*reaction 16: OH + OH <=> O + H2O */
    nuki[ 2 * kd + 15 ] += -1 ;
    nuki[ 2 * kd + 15 ] += -1 ;
    nuki[ 1 * kd + 15 ] += +1 ;
    nuki[ 6 * kd + 15 ] += +1 ;

    /*reaction 17: OH + H2 <=> H + H2O */
    nuki[ 2 * kd + 16 ] += -1 ;
    nuki[ 3 * kd + 16 ] += -1 ;
    nuki[ 0 * kd + 16 ] += +1 ;
    nuki[ 6 * kd + 16 ] += +1 ;

    /*reaction 18: H2 + O2 <=> HO2 + H */
    nuki[ 3 * kd + 17 ] += -1 ;
    nuki[ 4 * kd + 17 ] += -1 ;
    nuki[ 5 * kd + 17 ] += +1 ;
    nuki[ 0 * kd + 17 ] += +1 ;

    /*reaction 19: HO2 + H <=> OH + OH */
    nuki[ 5 * kd + 18 ] += -1 ;
    nuki[ 0 * kd + 18 ] += -1 ;
    nuki[ 2 * kd + 18 ] += +1 ;
    nuki[ 2 * kd + 18 ] += +1 ;

    /*reaction 20: HO2 + H <=> H2O + O */
    nuki[ 5 * kd + 19 ] += -1 ;
    nuki[ 0 * kd + 19 ] += -1 ;
    nuki[ 6 * kd + 19 ] += +1 ;
    nuki[ 1 * kd + 19 ] += +1 ;

    /*reaction 21: HO2 + O <=> OH + O2 */
    nuki[ 5 * kd + 20 ] += -1 ;
    nuki[ 1 * kd + 20 ] += -1 ;
    nuki[ 2 * kd + 20 ] += +1 ;
    nuki[ 4 * kd + 20 ] += +1 ;

    /*reaction 22: HO2 + OH <=> H2O + O2 */
    nuki[ 5 * kd + 21 ] += -1 ;
    nuki[ 2 * kd + 21 ] += -1 ;
    nuki[ 6 * kd + 21 ] += +1 ;
    nuki[ 4 * kd + 21 ] += +1 ;

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    nuki[ 5 * kd + 22 ] += -1 ;
    nuki[ 2 * kd + 22 ] += -1 ;
    nuki[ 6 * kd + 22 ] += +1 ;
    nuki[ 4 * kd + 22 ] += +1 ;

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    nuki[ 5 * kd + 23 ] += -1 ;
    nuki[ 2 * kd + 23 ] += -1 ;
    nuki[ 6 * kd + 23 ] += +1 ;
    nuki[ 4 * kd + 23 ] += +1 ;

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 5 * kd + 24 ] += -1 ;
    nuki[ 5 * kd + 24 ] += -1 ;
    nuki[ 7 * kd + 24 ] += +1 ;
    nuki[ 4 * kd + 24 ] += +1 ;

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 5 * kd + 25 ] += -1 ;
    nuki[ 5 * kd + 25 ] += -1 ;
    nuki[ 7 * kd + 25 ] += +1 ;
    nuki[ 4 * kd + 25 ] += +1 ;

    /*reaction 27: H2O2 + H <=> H2O + OH */
    nuki[ 7 * kd + 26 ] += -1 ;
    nuki[ 0 * kd + 26 ] += -1 ;
    nuki[ 6 * kd + 26 ] += +1 ;
    nuki[ 2 * kd + 26 ] += +1 ;

    /*reaction 28: H2O2 + H <=> HO2 + H2 */
    nuki[ 7 * kd + 27 ] += -1 ;
    nuki[ 0 * kd + 27 ] += -1 ;
    nuki[ 5 * kd + 27 ] += +1 ;
    nuki[ 3 * kd + 27 ] += +1 ;

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    nuki[ 7 * kd + 28 ] += -1 ;
    nuki[ 1 * kd + 28 ] += -1 ;
    nuki[ 5 * kd + 28 ] += +1 ;
    nuki[ 2 * kd + 28 ] += +1 ;

    /*reaction 30: H2O2 + OH <=> H2O + HO2 */
    nuki[ 7 * kd + 29 ] += -1 ;
    nuki[ 2 * kd + 29 ] += -1 ;
    nuki[ 6 * kd + 29 ] += +1 ;
    nuki[ 5 * kd + 29 ] += +1 ;

    /*reaction 31: H2O2 + OH <=> H2O + HO2 */
    nuki[ 7 * kd + 30 ] += -1 ;
    nuki[ 2 * kd + 30 ] += -1 ;
    nuki[ 6 * kd + 30 ] += +1 ;
    nuki[ 5 * kd + 30 ] += +1 ;

    /*reaction 32: NO + HO2 <=> NO2 + OH */
    nuki[ 8 * kd + 31 ] += -1 ;
    nuki[ 5 * kd + 31 ] += -1 ;
    nuki[ 9 * kd + 31 ] += +1 ;
    nuki[ 2 * kd + 31 ] += +1 ;

    /*reaction 33: NO2 + H <=> NO + OH */
    nuki[ 9 * kd + 32 ] += -1 ;
    nuki[ 0 * kd + 32 ] += -1 ;
    nuki[ 8 * kd + 32 ] += +1 ;
    nuki[ 2 * kd + 32 ] += +1 ;

    /*reaction 34: NO2 + O <=> NO + O2 */
    nuki[ 9 * kd + 33 ] += -1 ;
    nuki[ 1 * kd + 33 ] += -1 ;
    nuki[ 8 * kd + 33 ] += +1 ;
    nuki[ 4 * kd + 33 ] += +1 ;

    /*reaction 35: NO2 + NO2 <=> NO + NO + O2 */
    nuki[ 9 * kd + 34 ] += -1 ;
    nuki[ 9 * kd + 34 ] += -1 ;
    nuki[ 8 * kd + 34 ] += +1 ;
    nuki[ 8 * kd + 34 ] += +1 ;
    nuki[ 4 * kd + 34 ] += +1 ;

    /*reaction 36: N2O + H <=> N2 + OH */
    nuki[ 10 * kd + 35 ] += -1 ;
    nuki[ 0 * kd + 35 ] += -1 ;
    nuki[ 14 * kd + 35 ] += +1 ;
    nuki[ 2 * kd + 35 ] += +1 ;

    /*reaction 37: N2O + H <=> N2 + OH */
    nuki[ 10 * kd + 36 ] += -1 ;
    nuki[ 0 * kd + 36 ] += -1 ;
    nuki[ 14 * kd + 36 ] += +1 ;
    nuki[ 2 * kd + 36 ] += +1 ;

    /*reaction 38: N2O + O <=> NO + NO */
    nuki[ 10 * kd + 37 ] += -1 ;
    nuki[ 1 * kd + 37 ] += -1 ;
    nuki[ 8 * kd + 37 ] += +1 ;
    nuki[ 8 * kd + 37 ] += +1 ;

    /*reaction 39: N2O + O <=> N2 + O2 */
    nuki[ 10 * kd + 38 ] += -1 ;
    nuki[ 1 * kd + 38 ] += -1 ;
    nuki[ 14 * kd + 38 ] += +1 ;
    nuki[ 4 * kd + 38 ] += +1 ;

    /*reaction 40: NH + H <=> N + H2 */
    nuki[ 11 * kd + 39 ] += -1 ;
    nuki[ 0 * kd + 39 ] += -1 ;
    nuki[ 12 * kd + 39 ] += +1 ;
    nuki[ 3 * kd + 39 ] += +1 ;

    /*reaction 41: NH + O <=> NO + H */
    nuki[ 11 * kd + 40 ] += -1 ;
    nuki[ 1 * kd + 40 ] += -1 ;
    nuki[ 8 * kd + 40 ] += +1 ;
    nuki[ 0 * kd + 40 ] += +1 ;

    /*reaction 42: NH + OH <=> N + H2O */
    nuki[ 11 * kd + 41 ] += -1 ;
    nuki[ 2 * kd + 41 ] += -1 ;
    nuki[ 12 * kd + 41 ] += +1 ;
    nuki[ 6 * kd + 41 ] += +1 ;

    /*reaction 43: NH + O2 <=> NO + OH */
    nuki[ 11 * kd + 42 ] += -1 ;
    nuki[ 4 * kd + 42 ] += -1 ;
    nuki[ 8 * kd + 42 ] += +1 ;
    nuki[ 2 * kd + 42 ] += +1 ;

    /*reaction 44: NH + NO <=> N2O + H */
    nuki[ 11 * kd + 43 ] += -1 ;
    nuki[ 8 * kd + 43 ] += -1 ;
    nuki[ 10 * kd + 43 ] += +1 ;
    nuki[ 0 * kd + 43 ] += +1 ;

    /*reaction 45: NH + NO <=> N2O + H */
    nuki[ 11 * kd + 44 ] += -1 ;
    nuki[ 8 * kd + 44 ] += -1 ;
    nuki[ 10 * kd + 44 ] += +1 ;
    nuki[ 0 * kd + 44 ] += +1 ;

    /*reaction 46: NH + NO <=> N2 + OH */
    nuki[ 11 * kd + 45 ] += -1 ;
    nuki[ 8 * kd + 45 ] += -1 ;
    nuki[ 14 * kd + 45 ] += +1 ;
    nuki[ 2 * kd + 45 ] += +1 ;

    /*reaction 47: NH + NO2 <=> N2O + OH */
    nuki[ 11 * kd + 46 ] += -1 ;
    nuki[ 9 * kd + 46 ] += -1 ;
    nuki[ 10 * kd + 46 ] += +1 ;
    nuki[ 2 * kd + 46 ] += +1 ;

    /*reaction 48: N + OH <=> NO + H */
    nuki[ 12 * kd + 47 ] += -1 ;
    nuki[ 2 * kd + 47 ] += -1 ;
    nuki[ 8 * kd + 47 ] += +1 ;
    nuki[ 0 * kd + 47 ] += +1 ;

    /*reaction 49: N + O2 <=> NO + O */
    nuki[ 12 * kd + 48 ] += -1 ;
    nuki[ 4 * kd + 48 ] += -1 ;
    nuki[ 8 * kd + 48 ] += +1 ;
    nuki[ 1 * kd + 48 ] += +1 ;

    /*reaction 50: N + NO <=> N2 + O */
    nuki[ 12 * kd + 49 ] += -1 ;
    nuki[ 8 * kd + 49 ] += -1 ;
    nuki[ 14 * kd + 49 ] += +1 ;
    nuki[ 1 * kd + 49 ] += +1 ;

    /*reaction 51: NNH <=> N2 + H */
    nuki[ 13 * kd + 50 ] += -1 ;
    nuki[ 14 * kd + 50 ] += +1 ;
    nuki[ 0 * kd + 50 ] += +1 ;

    /*reaction 52: NNH + H <=> N2 + H2 */
    nuki[ 13 * kd + 51 ] += -1 ;
    nuki[ 0 * kd + 51 ] += -1 ;
    nuki[ 14 * kd + 51 ] += +1 ;
    nuki[ 3 * kd + 51 ] += +1 ;

    /*reaction 53: NNH + O <=> N2O + H */
    nuki[ 13 * kd + 52 ] += -1 ;
    nuki[ 1 * kd + 52 ] += -1 ;
    nuki[ 10 * kd + 52 ] += +1 ;
    nuki[ 0 * kd + 52 ] += +1 ;

    /*reaction 54: NNH + O <=> N2 + OH */
    nuki[ 13 * kd + 53 ] += -1 ;
    nuki[ 1 * kd + 53 ] += -1 ;
    nuki[ 14 * kd + 53 ] += +1 ;
    nuki[ 2 * kd + 53 ] += +1 ;

    /*reaction 55: NNH + O <=> NH + NO */
    nuki[ 13 * kd + 54 ] += -1 ;
    nuki[ 1 * kd + 54 ] += -1 ;
    nuki[ 11 * kd + 54 ] += +1 ;
    nuki[ 8 * kd + 54 ] += +1 ;

    /*reaction 56: NNH + OH <=> N2 + H2O */
    nuki[ 13 * kd + 55 ] += -1 ;
    nuki[ 2 * kd + 55 ] += -1 ;
    nuki[ 14 * kd + 55 ] += +1 ;
    nuki[ 6 * kd + 55 ] += +1 ;

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    nuki[ 13 * kd + 56 ] += -1 ;
    nuki[ 4 * kd + 56 ] += -1 ;
    nuki[ 14 * kd + 56 ] += +1 ;
    nuki[ 5 * kd + 56 ] += +1 ;

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    nuki[ 13 * kd + 57 ] += -1 ;
    nuki[ 4 * kd + 57 ] += -1 ;
    nuki[ 14 * kd + 57 ] += +1 ;
    nuki[ 0 * kd + 57 ] += +1 ;
    nuki[ 4 * kd + 57 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 15; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 1 * kd + 0 ] = 1; /*O */

    /*OH */
    ncf[ 2 * kd + 0 ] = 1; /*O */
    ncf[ 2 * kd + 1 ] = 1; /*H */

    /*H2 */
    ncf[ 3 * kd + 1 ] = 2; /*H */

    /*O2 */
    ncf[ 4 * kd + 0 ] = 2; /*O */

    /*HO2 */
    ncf[ 5 * kd + 1 ] = 1; /*H */
    ncf[ 5 * kd + 0 ] = 2; /*O */

    /*H2O */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 0 ] = 1; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*NO */
    ncf[ 8 * kd + 2 ] = 1; /*N */
    ncf[ 8 * kd + 0 ] = 1; /*O */

    /*NO2 */
    ncf[ 9 * kd + 2 ] = 1; /*N */
    ncf[ 9 * kd + 0 ] = 2; /*O */

    /*N2O */
    ncf[ 10 * kd + 2 ] = 2; /*N */
    ncf[ 10 * kd + 0 ] = 1; /*O */

    /*NH */
    ncf[ 11 * kd + 2 ] = 1; /*N */
    ncf[ 11 * kd + 1 ] = 1; /*H */

    /*N */
    ncf[ 12 * kd + 2 ] = 1; /*N */

    /*NNH */
    ncf[ 13 * kd + 2 ] = 2; /*N */
    ncf[ 13 * kd + 1 ] = 1; /*H */

    /*N2 */
    ncf[ 14 * kd + 2 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)
{
    for (int i=0; i<58; ++i) {
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
    double gort[15]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + O2 (+N2) <=> HO2 (+N2) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[2] *= 1e-06; 

    /*reaction 4: NO + O (+M) <=> NO2 (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: N2O (+M) <=> N2 + O (+M) */
    eqcon[4] *= 1e-06; 

    /*reaction 6: H + H + M <=> H2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O + M <=> OH + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: O + O + M <=> O2 + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: OH + H + M <=> H2O + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 <=> O + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H + H + N2 <=> H2 + N2 */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + H + H2 <=> H2 + H2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + H + H2O <=> H2 + H2O */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + H2 <=> OH + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + H2 <=> OH + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: OH + OH <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + H2 <=> H + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2 + O2 <=> HO2 + H */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + H <=> H2O + O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + O <=> OH + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + OH <=> H2O + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + H <=> H2O + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + H <=> HO2 + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: NO + HO2 <=> NO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: NO2 + H <=> NO + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: NO2 + O <=> NO + O2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: NO2 + NO2 <=> NO + NO + O2 */
    eqcon[34] *= 1e-06; 

    /*reaction 36: N2O + H <=> N2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: N2O + O <=> NO + NO */
    /*eqcon[37] *= 1;  */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + H <=> N + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + O <=> NO + H */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> N + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + NO <=> N2O + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + NO <=> N2O + H */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*eqcon[46] *= 1;  */

    /*reaction 48: N + OH <=> NO + H */
    /*eqcon[47] *= 1;  */

    /*reaction 49: N + O2 <=> NO + O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: N + NO <=> N2 + O */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NNH <=> N2 + H */
    eqcon[50] *= 1e-06; 

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NNH + O <=> N2O + H */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*eqcon[53] *= 1;  */

    /*reaction 55: NNH + O <=> NH + NO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    eqcon[57] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[15]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + O2 (+N2) <=> HO2 (+N2) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[2] *= 1e-06; 

    /*reaction 4: NO + O (+M) <=> NO2 (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: N2O (+M) <=> N2 + O (+M) */
    eqcon[4] *= 1e-06; 

    /*reaction 6: H + H + M <=> H2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O + M <=> OH + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: O + O + M <=> O2 + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: OH + H + M <=> H2O + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 <=> O + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H + H + N2 <=> H2 + N2 */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + H + H2 <=> H2 + H2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + H + H2O <=> H2 + H2O */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + H2 <=> OH + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + H2 <=> OH + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: OH + OH <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + H2 <=> H + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2 + O2 <=> HO2 + H */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + H <=> H2O + O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + O <=> OH + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + OH <=> H2O + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + H <=> H2O + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + H <=> HO2 + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: NO + HO2 <=> NO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: NO2 + H <=> NO + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: NO2 + O <=> NO + O2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: NO2 + NO2 <=> NO + NO + O2 */
    eqcon[34] *= 1e-06; 

    /*reaction 36: N2O + H <=> N2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: N2O + O <=> NO + NO */
    /*eqcon[37] *= 1;  */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + H <=> N + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + O <=> NO + H */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> N + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + NO <=> N2O + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + NO <=> N2O + H */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*eqcon[46] *= 1;  */

    /*reaction 48: N + OH <=> NO + H */
    /*eqcon[47] *= 1;  */

    /*reaction 49: N + O2 <=> NO + O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: N + NO <=> N2 + O */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NNH <=> N2 + H */
    eqcon[50] *= 1e-06; 

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NNH + O <=> N2O + H */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*eqcon[53] *= 1;  */

    /*reaction 55: NNH + O <=> NH + NO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    eqcon[57] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[15]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + O2 (+N2) <=> HO2 (+N2) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[2] *= 1e-06; 

    /*reaction 4: NO + O (+M) <=> NO2 (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: N2O (+M) <=> N2 + O (+M) */
    eqcon[4] *= 1e-06; 

    /*reaction 6: H + H + M <=> H2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O + M <=> OH + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: O + O + M <=> O2 + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: OH + H + M <=> H2O + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 <=> O + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H + H + N2 <=> H2 + N2 */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + H + H2 <=> H2 + H2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + H + H2O <=> H2 + H2O */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + H2 <=> OH + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + H2 <=> OH + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: OH + OH <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + H2 <=> H + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2 + O2 <=> HO2 + H */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + H <=> H2O + O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + O <=> OH + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + OH <=> H2O + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + H <=> H2O + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + H <=> HO2 + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: NO + HO2 <=> NO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: NO2 + H <=> NO + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: NO2 + O <=> NO + O2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: NO2 + NO2 <=> NO + NO + O2 */
    eqcon[34] *= 1e-06; 

    /*reaction 36: N2O + H <=> N2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: N2O + O <=> NO + NO */
    /*eqcon[37] *= 1;  */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + H <=> N + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + O <=> NO + H */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> N + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + NO <=> N2O + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + NO <=> N2O + H */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*eqcon[46] *= 1;  */

    /*reaction 48: N + OH <=> NO + H */
    /*eqcon[47] *= 1;  */

    /*reaction 49: N + O2 <=> NO + O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: N + NO <=> N2 + O */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NNH <=> N2 + H */
    eqcon[50] *= 1e-06; 

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NNH + O <=> N2O + H */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*eqcon[53] *= 1;  */

    /*reaction 55: NNH + O <=> NH + NO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    eqcon[57] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[15]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + O2 (+N2) <=> HO2 (+N2) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[2] *= 1e-06; 

    /*reaction 4: NO + O (+M) <=> NO2 (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: N2O (+M) <=> N2 + O (+M) */
    eqcon[4] *= 1e-06; 

    /*reaction 6: H + H + M <=> H2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O + M <=> OH + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: O + O + M <=> O2 + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: OH + H + M <=> H2O + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 <=> O + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H + H + N2 <=> H2 + N2 */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + H + H2 <=> H2 + H2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + H + H2O <=> H2 + H2O */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + H2 <=> OH + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + H2 <=> OH + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: OH + OH <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + H2 <=> H + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2 + O2 <=> HO2 + H */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + H <=> H2O + O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + O <=> OH + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + OH <=> H2O + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + H <=> H2O + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + H <=> HO2 + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: NO + HO2 <=> NO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: NO2 + H <=> NO + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: NO2 + O <=> NO + O2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: NO2 + NO2 <=> NO + NO + O2 */
    eqcon[34] *= 1e-06; 

    /*reaction 36: N2O + H <=> N2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: N2O + O <=> NO + NO */
    /*eqcon[37] *= 1;  */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + H <=> N + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + O <=> NO + H */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> N + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + NO <=> N2O + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + NO <=> N2O + H */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*eqcon[46] *= 1;  */

    /*reaction 48: N + OH <=> NO + H */
    /*eqcon[47] *= 1;  */

    /*reaction 49: N + O2 <=> NO + O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: N + NO <=> N2 + O */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NNH <=> N2 + H */
    eqcon[50] *= 1e-06; 

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NNH + O <=> N2O + H */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*eqcon[53] *= 1;  */

    /*reaction 55: NNH + O <=> NH + NO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    eqcon[57] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[15]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + O2 (+N2) <=> HO2 (+N2) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[2] *= 1e-06; 

    /*reaction 4: NO + O (+M) <=> NO2 (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: N2O (+M) <=> N2 + O (+M) */
    eqcon[4] *= 1e-06; 

    /*reaction 6: H + H + M <=> H2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O + M <=> OH + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: O + O + M <=> O2 + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: OH + H + M <=> H2O + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 <=> O + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H + H + N2 <=> H2 + N2 */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + H + H2 <=> H2 + H2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + H + H2O <=> H2 + H2O */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + H2 <=> OH + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + H2 <=> OH + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: OH + OH <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + H2 <=> H + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2 + O2 <=> HO2 + H */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + H <=> H2O + O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + O <=> OH + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + OH <=> H2O + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + H <=> H2O + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + H <=> HO2 + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: NO + HO2 <=> NO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: NO2 + H <=> NO + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: NO2 + O <=> NO + O2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: NO2 + NO2 <=> NO + NO + O2 */
    eqcon[34] *= 1e-06; 

    /*reaction 36: N2O + H <=> N2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: N2O + O <=> NO + NO */
    /*eqcon[37] *= 1;  */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + H <=> N + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + O <=> NO + H */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> N + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + NO <=> N2O + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + NO <=> N2O + H */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*eqcon[46] *= 1;  */

    /*reaction 48: N + OH <=> NO + H */
    /*eqcon[47] *= 1;  */

    /*reaction 49: N + O2 <=> NO + O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: N + NO <=> N2 + O */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NNH <=> N2 + H */
    eqcon[50] *= 1e-06; 

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NNH + O <=> N2O + H */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*eqcon[53] *= 1;  */

    /*reaction 55: NNH + O <=> NH + NO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    eqcon[57] *= 1e-06; 
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[58];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[58];
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

    double qdot, q_f[58], q_r[58];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 15; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[0] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[0] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[2] += qdot;
    wdot[2] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] += qdot;
    wdot[10] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[0] -= qdot;
    wdot[0] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[1] -= qdot;
    wdot[2] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[1] -= qdot;
    wdot[1] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] += qdot;
    wdot[4] -= qdot;

    qdot = q_f[10]-q_r[10];
    wdot[0] -= qdot;
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[14] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[0] -= qdot;
    wdot[0] -= qdot;
    wdot[3] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[0] -= qdot;
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;

    qdot = q_f[14]-q_r[14];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;

    qdot = q_f[15]-q_r[15];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[0] += qdot;
    wdot[3] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[2] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[20]-q_r[20];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[21]-q_r[21];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[23]-q_r[23];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[24]-q_r[24];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[5] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[25]-q_r[25];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[5] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[26]-q_r[26];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[29]-q_r[29];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[30]-q_r[30];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[31]-q_r[31];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[8] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[33]-q_r[33];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[8] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[34]-q_r[34];
    wdot[4] += qdot;
    wdot[8] += qdot;
    wdot[8] += qdot;
    wdot[9] -= qdot;
    wdot[9] -= qdot;

    qdot = q_f[35]-q_r[35];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[10] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[10] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[37]-q_r[37];
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[38]-q_r[38];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[10] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[39]-q_r[39];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[40]-q_r[40];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[41]-q_r[41];
    wdot[2] -= qdot;
    wdot[6] += qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[42]-q_r[42];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[43]-q_r[43];
    wdot[0] += qdot;
    wdot[8] -= qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[44]-q_r[44];
    wdot[0] += qdot;
    wdot[8] -= qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[45]-q_r[45];
    wdot[2] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[46]-q_r[46];
    wdot[2] += qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[47]-q_r[47];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[8] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[48]-q_r[48];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[8] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[12] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[50]-q_r[50];
    wdot[0] += qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[51]-q_r[51];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[52]-q_r[52];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[10] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[53]-q_r[53];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[54]-q_r[54];
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[55]-q_r[55];
    wdot[2] -= qdot;
    wdot[6] += qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[56]-q_r[56];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[57]-q_r[57];
    wdot[0] += qdot;
    wdot[4] -= qdot;
    wdot[4] += qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    return;
}

void comp_k_f(double * restrict tc, double invT, double * restrict k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<58; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double * restrict tc, double invT, double * restrict Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[15];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[0] + g_RT[4] - g_RT[5];
    Kc[1] = g_RT[0] + g_RT[4] - g_RT[5];
    Kc[2] = -g_RT[2] - g_RT[2] + g_RT[7];
    Kc[3] = g_RT[1] + g_RT[8] - g_RT[9];
    Kc[4] = -g_RT[1] + g_RT[10] - g_RT[14];
    Kc[5] = g_RT[0] + g_RT[0] - g_RT[3];
    Kc[6] = g_RT[0] + g_RT[1] - g_RT[2];
    Kc[7] = g_RT[1] + g_RT[1] - g_RT[4];
    Kc[8] = g_RT[0] + g_RT[2] - g_RT[6];
    Kc[9] = g_RT[0] - g_RT[1] - g_RT[2] + g_RT[4];
    Kc[10] = g_RT[0] + g_RT[0] - g_RT[3] + g_RT[14] - g_RT[14];
    Kc[11] = g_RT[0] + g_RT[0] + g_RT[3] - g_RT[3] - g_RT[3];
    Kc[12] = g_RT[0] + g_RT[0] - g_RT[3] + g_RT[6] - g_RT[6];
    Kc[13] = -g_RT[0] + g_RT[1] - g_RT[2] + g_RT[3];
    Kc[14] = -g_RT[0] + g_RT[1] - g_RT[2] + g_RT[3];
    Kc[15] = -g_RT[1] + g_RT[2] + g_RT[2] - g_RT[6];
    Kc[16] = -g_RT[0] + g_RT[2] + g_RT[3] - g_RT[6];
    Kc[17] = -g_RT[0] + g_RT[3] + g_RT[4] - g_RT[5];
    Kc[18] = g_RT[0] - g_RT[2] - g_RT[2] + g_RT[5];
    Kc[19] = g_RT[0] - g_RT[1] + g_RT[5] - g_RT[6];
    Kc[20] = g_RT[1] - g_RT[2] - g_RT[4] + g_RT[5];
    Kc[21] = g_RT[2] - g_RT[4] + g_RT[5] - g_RT[6];
    Kc[22] = g_RT[2] - g_RT[4] + g_RT[5] - g_RT[6];
    Kc[23] = g_RT[2] - g_RT[4] + g_RT[5] - g_RT[6];
    Kc[24] = -g_RT[4] + g_RT[5] + g_RT[5] - g_RT[7];
    Kc[25] = -g_RT[4] + g_RT[5] + g_RT[5] - g_RT[7];
    Kc[26] = g_RT[0] - g_RT[2] - g_RT[6] + g_RT[7];
    Kc[27] = g_RT[0] - g_RT[3] - g_RT[5] + g_RT[7];
    Kc[28] = g_RT[1] - g_RT[2] - g_RT[5] + g_RT[7];
    Kc[29] = g_RT[2] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[30] = g_RT[2] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[31] = -g_RT[2] + g_RT[5] + g_RT[8] - g_RT[9];
    Kc[32] = g_RT[0] - g_RT[2] - g_RT[8] + g_RT[9];
    Kc[33] = g_RT[1] - g_RT[4] - g_RT[8] + g_RT[9];
    Kc[34] = -g_RT[4] - g_RT[8] - g_RT[8] + g_RT[9] + g_RT[9];
    Kc[35] = g_RT[0] - g_RT[2] + g_RT[10] - g_RT[14];
    Kc[36] = g_RT[0] - g_RT[2] + g_RT[10] - g_RT[14];
    Kc[37] = g_RT[1] - g_RT[8] - g_RT[8] + g_RT[10];
    Kc[38] = g_RT[1] - g_RT[4] + g_RT[10] - g_RT[14];
    Kc[39] = g_RT[0] - g_RT[3] + g_RT[11] - g_RT[12];
    Kc[40] = -g_RT[0] + g_RT[1] - g_RT[8] + g_RT[11];
    Kc[41] = g_RT[2] - g_RT[6] + g_RT[11] - g_RT[12];
    Kc[42] = -g_RT[2] + g_RT[4] - g_RT[8] + g_RT[11];
    Kc[43] = -g_RT[0] + g_RT[8] - g_RT[10] + g_RT[11];
    Kc[44] = -g_RT[0] + g_RT[8] - g_RT[10] + g_RT[11];
    Kc[45] = -g_RT[2] + g_RT[8] + g_RT[11] - g_RT[14];
    Kc[46] = -g_RT[2] + g_RT[9] - g_RT[10] + g_RT[11];
    Kc[47] = -g_RT[0] + g_RT[2] - g_RT[8] + g_RT[12];
    Kc[48] = -g_RT[1] + g_RT[4] - g_RT[8] + g_RT[12];
    Kc[49] = -g_RT[1] + g_RT[8] + g_RT[12] - g_RT[14];
    Kc[50] = -g_RT[0] + g_RT[13] - g_RT[14];
    Kc[51] = g_RT[0] - g_RT[3] + g_RT[13] - g_RT[14];
    Kc[52] = -g_RT[0] + g_RT[1] - g_RT[10] + g_RT[13];
    Kc[53] = g_RT[1] - g_RT[2] + g_RT[13] - g_RT[14];
    Kc[54] = g_RT[1] - g_RT[8] - g_RT[11] + g_RT[13];
    Kc[55] = g_RT[2] - g_RT[6] + g_RT[13] - g_RT[14];
    Kc[56] = g_RT[4] - g_RT[5] + g_RT[13] - g_RT[14];
    Kc[57] = -g_RT[0] + g_RT[4] - g_RT[4] + g_RT[13] - g_RT[14];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<58; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refCinv;
    Kc[2] *= refC;
    Kc[3] *= refCinv;
    Kc[4] *= refC;
    Kc[5] *= refCinv;
    Kc[6] *= refCinv;
    Kc[7] *= refCinv;
    Kc[8] *= refCinv;
    Kc[10] *= refCinv;
    Kc[11] *= refCinv;
    Kc[12] *= refCinv;
    Kc[34] *= refC;
    Kc[50] *= refC;
    Kc[57] *= refC;

    return;
}

void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[0]*sc[4];
    qr[0] = sc[5];

    /*reaction 2: H + O2 (+N2) <=> HO2 (+N2) */
    qf[1] = sc[0]*sc[4];
    qr[1] = sc[5];

    /*reaction 3: H2O2 (+M) <=> OH + OH (+M) */
    qf[2] = sc[7];
    qr[2] = sc[2]*sc[2];

    /*reaction 4: NO + O (+M) <=> NO2 (+M) */
    qf[3] = sc[1]*sc[8];
    qr[3] = sc[9];

    /*reaction 5: N2O (+M) <=> N2 + O (+M) */
    qf[4] = sc[10];
    qr[4] = sc[1]*sc[14];

    /*reaction 6: H + H + M <=> H2 + M */
    qf[5] = sc[0]*sc[0];
    qr[5] = sc[3];

    /*reaction 7: H + O + M <=> OH + M */
    qf[6] = sc[0]*sc[1];
    qr[6] = sc[2];

    /*reaction 8: O + O + M <=> O2 + M */
    qf[7] = sc[1]*sc[1];
    qr[7] = sc[4];

    /*reaction 9: OH + H + M <=> H2O + M */
    qf[8] = sc[0]*sc[2];
    qr[8] = sc[6];

    /*reaction 10: H + O2 <=> O + OH */
    qf[9] = sc[0]*sc[4];
    qr[9] = sc[1]*sc[2];

    /*reaction 11: H + H + N2 <=> H2 + N2 */
    qf[10] = sc[0]*sc[0]*sc[14];
    qr[10] = sc[3]*sc[14];

    /*reaction 12: H + H + H2 <=> H2 + H2 */
    qf[11] = sc[0]*sc[0]*sc[3];
    qr[11] = sc[3]*sc[3];

    /*reaction 13: H + H + H2O <=> H2 + H2O */
    qf[12] = sc[0]*sc[0]*sc[6];
    qr[12] = sc[3]*sc[6];

    /*reaction 14: O + H2 <=> OH + H */
    qf[13] = sc[1]*sc[3];
    qr[13] = sc[0]*sc[2];

    /*reaction 15: O + H2 <=> OH + H */
    qf[14] = sc[1]*sc[3];
    qr[14] = sc[0]*sc[2];

    /*reaction 16: OH + OH <=> O + H2O */
    qf[15] = sc[2]*sc[2];
    qr[15] = sc[1]*sc[6];

    /*reaction 17: OH + H2 <=> H + H2O */
    qf[16] = sc[2]*sc[3];
    qr[16] = sc[0]*sc[6];

    /*reaction 18: H2 + O2 <=> HO2 + H */
    qf[17] = sc[3]*sc[4];
    qr[17] = sc[0]*sc[5];

    /*reaction 19: HO2 + H <=> OH + OH */
    qf[18] = sc[0]*sc[5];
    qr[18] = sc[2]*sc[2];

    /*reaction 20: HO2 + H <=> H2O + O */
    qf[19] = sc[0]*sc[5];
    qr[19] = sc[1]*sc[6];

    /*reaction 21: HO2 + O <=> OH + O2 */
    qf[20] = sc[1]*sc[5];
    qr[20] = sc[2]*sc[4];

    /*reaction 22: HO2 + OH <=> H2O + O2 */
    qf[21] = sc[2]*sc[5];
    qr[21] = sc[4]*sc[6];

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    qf[22] = sc[2]*sc[5];
    qr[22] = sc[4]*sc[6];

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    qf[23] = sc[2]*sc[5];
    qr[23] = sc[4]*sc[6];

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    qf[24] = sc[5]*sc[5];
    qr[24] = sc[4]*sc[7];

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    qf[25] = sc[5]*sc[5];
    qr[25] = sc[4]*sc[7];

    /*reaction 27: H2O2 + H <=> H2O + OH */
    qf[26] = sc[0]*sc[7];
    qr[26] = sc[2]*sc[6];

    /*reaction 28: H2O2 + H <=> HO2 + H2 */
    qf[27] = sc[0]*sc[7];
    qr[27] = sc[3]*sc[5];

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    qf[28] = sc[1]*sc[7];
    qr[28] = sc[2]*sc[5];

    /*reaction 30: H2O2 + OH <=> H2O + HO2 */
    qf[29] = sc[2]*sc[7];
    qr[29] = sc[5]*sc[6];

    /*reaction 31: H2O2 + OH <=> H2O + HO2 */
    qf[30] = sc[2]*sc[7];
    qr[30] = sc[5]*sc[6];

    /*reaction 32: NO + HO2 <=> NO2 + OH */
    qf[31] = sc[5]*sc[8];
    qr[31] = sc[2]*sc[9];

    /*reaction 33: NO2 + H <=> NO + OH */
    qf[32] = sc[0]*sc[9];
    qr[32] = sc[2]*sc[8];

    /*reaction 34: NO2 + O <=> NO + O2 */
    qf[33] = sc[1]*sc[9];
    qr[33] = sc[4]*sc[8];

    /*reaction 35: NO2 + NO2 <=> NO + NO + O2 */
    qf[34] = sc[9]*sc[9];
    qr[34] = sc[4]*sc[8]*sc[8];

    /*reaction 36: N2O + H <=> N2 + OH */
    qf[35] = sc[0]*sc[10];
    qr[35] = sc[2]*sc[14];

    /*reaction 37: N2O + H <=> N2 + OH */
    qf[36] = sc[0]*sc[10];
    qr[36] = sc[2]*sc[14];

    /*reaction 38: N2O + O <=> NO + NO */
    qf[37] = sc[1]*sc[10];
    qr[37] = sc[8]*sc[8];

    /*reaction 39: N2O + O <=> N2 + O2 */
    qf[38] = sc[1]*sc[10];
    qr[38] = sc[4]*sc[14];

    /*reaction 40: NH + H <=> N + H2 */
    qf[39] = sc[0]*sc[11];
    qr[39] = sc[3]*sc[12];

    /*reaction 41: NH + O <=> NO + H */
    qf[40] = sc[1]*sc[11];
    qr[40] = sc[0]*sc[8];

    /*reaction 42: NH + OH <=> N + H2O */
    qf[41] = sc[2]*sc[11];
    qr[41] = sc[6]*sc[12];

    /*reaction 43: NH + O2 <=> NO + OH */
    qf[42] = sc[4]*sc[11];
    qr[42] = sc[2]*sc[8];

    /*reaction 44: NH + NO <=> N2O + H */
    qf[43] = sc[8]*sc[11];
    qr[43] = sc[0]*sc[10];

    /*reaction 45: NH + NO <=> N2O + H */
    qf[44] = sc[8]*sc[11];
    qr[44] = sc[0]*sc[10];

    /*reaction 46: NH + NO <=> N2 + OH */
    qf[45] = sc[8]*sc[11];
    qr[45] = sc[2]*sc[14];

    /*reaction 47: NH + NO2 <=> N2O + OH */
    qf[46] = sc[9]*sc[11];
    qr[46] = sc[2]*sc[10];

    /*reaction 48: N + OH <=> NO + H */
    qf[47] = sc[2]*sc[12];
    qr[47] = sc[0]*sc[8];

    /*reaction 49: N + O2 <=> NO + O */
    qf[48] = sc[4]*sc[12];
    qr[48] = sc[1]*sc[8];

    /*reaction 50: N + NO <=> N2 + O */
    qf[49] = sc[8]*sc[12];
    qr[49] = sc[1]*sc[14];

    /*reaction 51: NNH <=> N2 + H */
    qf[50] = sc[13];
    qr[50] = sc[0]*sc[14];

    /*reaction 52: NNH + H <=> N2 + H2 */
    qf[51] = sc[0]*sc[13];
    qr[51] = sc[3]*sc[14];

    /*reaction 53: NNH + O <=> N2O + H */
    qf[52] = sc[1]*sc[13];
    qr[52] = sc[0]*sc[10];

    /*reaction 54: NNH + O <=> N2 + OH */
    qf[53] = sc[1]*sc[13];
    qr[53] = sc[2]*sc[14];

    /*reaction 55: NNH + O <=> NH + NO */
    qf[54] = sc[1]*sc[13];
    qr[54] = sc[8]*sc[11];

    /*reaction 56: NNH + OH <=> N2 + H2O */
    qf[55] = sc[2]*sc[13];
    qr[55] = sc[6]*sc[14];

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    qf[56] = sc[4]*sc[13];
    qr[56] = sc[5]*sc[14];

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    qf[57] = sc[4]*sc[13];
    qr[57] = sc[0]*sc[4]*sc[14];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 15; ++i) {
        mixture += sc[i];
    }

    double Corr[58];
    for (int i = 0; i < 58; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[4];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[14] + (TB[0][1] - 1)*sc[6] + (TB[0][2] - 1)*sc[3] + (TB[0][3] - 1)*sc[4];
        alpha[1] = sc[14];
        alpha[2] = mixture + (TB[2][0] - 1)*sc[6] + (TB[2][1] - 1)*sc[3];
        alpha[3] = mixture;
#ifdef __INTEL_COMPILER
         #pragma simd
#endif
        for (int i=0; i<4; i++)
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
        alpha = mixture + (TB[4][0] - 1)*sc[14] + (TB[4][1] - 1)*sc[4] + (TB[4][2] - 1)*sc[6];
        double redP = alpha / k_f_save[4] * phase_units[4] * low_A[4] * exp(low_beta[4] * tc[0] - activation_units[4] * low_Ea[4] * invT);
        Corr[4] = redP / (1. + redP);
    }

    /* simple three-body correction */
    {
        double alpha;
        alpha = mixture + (TB[5][0] - 1)*sc[14] + (TB[5][1] - 1)*sc[6] + (TB[5][2] - 1)*sc[3];
        Corr[5] = alpha;
        alpha = mixture + (TB[6][0] - 1)*sc[6];
        Corr[6] = alpha;
        alpha = mixture + (TB[7][0] - 1)*sc[14] + (TB[7][1] - 1)*sc[4] + (TB[7][2] - 1)*sc[6];
        Corr[7] = alpha;
        alpha = mixture + (TB[8][0] - 1)*sc[3] + (TB[8][1] - 1)*sc[6];
        Corr[8] = alpha;
    }

    for (int i=0; i<58; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[58*npt], Kc_s[58*npt], mixture[npt], g_RT[15*npt];
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

    for (int n=0; n<15; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

    vcomp_k_f(npt, k_f_s, tc, invT);

    vcomp_gibbs(npt, g_RT, tc);

    vcomp_Kc(npt, Kc_s, g_RT, invT);

    vcomp_wdot_1_50(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
    vcomp_wdot_51_58(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
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
    }
}

void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[15];
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

        Kc_s[0*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[2*npt+i] = refC * exp((g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[2*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[8*npt+i]) - (g_RT[9*npt+i]));
        Kc_s[4*npt+i] = refC * exp((g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[0*npt+i]) - (g_RT[3*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[1*npt+i]) - (g_RT[2*npt+i]));
        Kc_s[7*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[1*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[8*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[9*npt+i] = exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i] + g_RT[2*npt+i]));
        Kc_s[10*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[0*npt+i] + g_RT[14*npt+i]) - (g_RT[3*npt+i] + g_RT[14*npt+i]));
        Kc_s[11*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[0*npt+i] + g_RT[3*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[12*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[0*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[13*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[0*npt+i] + g_RT[2*npt+i]));
        Kc_s[14*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[0*npt+i] + g_RT[2*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[2*npt+i] + g_RT[2*npt+i]) - (g_RT[1*npt+i] + g_RT[6*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[2*npt+i] + g_RT[3*npt+i]) - (g_RT[0*npt+i] + g_RT[6*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[3*npt+i] + g_RT[4*npt+i]) - (g_RT[0*npt+i] + g_RT[5*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[2*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[1*npt+i] + g_RT[6*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[2*npt+i] + g_RT[5*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[2*npt+i] + g_RT[5*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[2*npt+i] + g_RT[5*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[5*npt+i] + g_RT[5*npt+i]) - (g_RT[4*npt+i] + g_RT[7*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[5*npt+i] + g_RT[5*npt+i]) - (g_RT[4*npt+i] + g_RT[7*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[0*npt+i] + g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[6*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[0*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[29*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[31*npt+i] = exp((g_RT[5*npt+i] + g_RT[8*npt+i]) - (g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[0*npt+i] + g_RT[9*npt+i]) - (g_RT[2*npt+i] + g_RT[8*npt+i]));
        Kc_s[33*npt+i] = exp((g_RT[1*npt+i] + g_RT[9*npt+i]) - (g_RT[4*npt+i] + g_RT[8*npt+i]));
        Kc_s[34*npt+i] = refC * exp((g_RT[9*npt+i] + g_RT[9*npt+i]) - (g_RT[4*npt+i] + g_RT[8*npt+i] + g_RT[8*npt+i]));
        Kc_s[35*npt+i] = exp((g_RT[0*npt+i] + g_RT[10*npt+i]) - (g_RT[2*npt+i] + g_RT[14*npt+i]));
        Kc_s[36*npt+i] = exp((g_RT[0*npt+i] + g_RT[10*npt+i]) - (g_RT[2*npt+i] + g_RT[14*npt+i]));
        Kc_s[37*npt+i] = exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[8*npt+i] + g_RT[8*npt+i]));
        Kc_s[38*npt+i] = exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[4*npt+i] + g_RT[14*npt+i]));
        Kc_s[39*npt+i] = exp((g_RT[0*npt+i] + g_RT[11*npt+i]) - (g_RT[3*npt+i] + g_RT[12*npt+i]));
        Kc_s[40*npt+i] = exp((g_RT[1*npt+i] + g_RT[11*npt+i]) - (g_RT[0*npt+i] + g_RT[8*npt+i]));
        Kc_s[41*npt+i] = exp((g_RT[2*npt+i] + g_RT[11*npt+i]) - (g_RT[6*npt+i] + g_RT[12*npt+i]));
        Kc_s[42*npt+i] = exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[8*npt+i]));
        Kc_s[43*npt+i] = exp((g_RT[8*npt+i] + g_RT[11*npt+i]) - (g_RT[0*npt+i] + g_RT[10*npt+i]));
        Kc_s[44*npt+i] = exp((g_RT[8*npt+i] + g_RT[11*npt+i]) - (g_RT[0*npt+i] + g_RT[10*npt+i]));
        Kc_s[45*npt+i] = exp((g_RT[8*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[14*npt+i]));
        Kc_s[46*npt+i] = exp((g_RT[9*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[10*npt+i]));
        Kc_s[47*npt+i] = exp((g_RT[2*npt+i] + g_RT[12*npt+i]) - (g_RT[0*npt+i] + g_RT[8*npt+i]));
        Kc_s[48*npt+i] = exp((g_RT[4*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[8*npt+i]));
        Kc_s[49*npt+i] = exp((g_RT[8*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[50*npt+i] = refC * exp((g_RT[13*npt+i]) - (g_RT[0*npt+i] + g_RT[14*npt+i]));
        Kc_s[51*npt+i] = exp((g_RT[0*npt+i] + g_RT[13*npt+i]) - (g_RT[3*npt+i] + g_RT[14*npt+i]));
        Kc_s[52*npt+i] = exp((g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[0*npt+i] + g_RT[10*npt+i]));
        Kc_s[53*npt+i] = exp((g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[2*npt+i] + g_RT[14*npt+i]));
        Kc_s[54*npt+i] = exp((g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[8*npt+i] + g_RT[11*npt+i]));
        Kc_s[55*npt+i] = exp((g_RT[2*npt+i] + g_RT[13*npt+i]) - (g_RT[6*npt+i] + g_RT[14*npt+i]));
        Kc_s[56*npt+i] = exp((g_RT[4*npt+i] + g_RT[13*npt+i]) - (g_RT[5*npt+i] + g_RT[14*npt+i]));
        Kc_s[57*npt+i] = refC * exp((g_RT[4*npt+i] + g_RT[13*npt+i]) - (g_RT[0*npt+i] + g_RT[4*npt+i] + g_RT[14*npt+i]));
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

        /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[14*npt+i] + (TB[0][1] - 1)*sc[6*npt+i] + (TB[0][2] - 1)*sc[3*npt+i] + (TB[0][3] - 1)*sc[4*npt+i];
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
        phi_r = sc[5*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 2: H + O2 (+N2) <=> HO2 (+N2) */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        alpha = sc[14*npt+i];
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
        phi_r = sc[5*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 3: H2O2 (+M) <=> OH + OH (+M) */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[6*npt+i] + (TB[2][1] - 1)*sc[3*npt+i];
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
        phi_r = sc[2*npt+i]*sc[2*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[2*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 4: NO + O (+M) <=> NO2 (+M) */
        phi_f = sc[1*npt+i]*sc[8*npt+i];
        alpha = mixture[i];
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
        phi_r = sc[9*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 5: N2O (+M) <=> N2 + O (+M) */
        phi_f = sc[10*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[14*npt+i] + (TB[4][1] - 1)*sc[4*npt+i] + (TB[4][2] - 1)*sc[6*npt+i];
        k_f = k_f_s[4*npt+i];
        redP = alpha / k_f * phase_units[4] * low_A[4] * exp(low_beta[4] * tc[i] - activation_units[4] * low_Ea[4] * invT[i]);
        F = redP / (1 + redP);
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 6: H + H + M <=> H2 + M */
        phi_f = sc[0*npt+i]*sc[0*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[14*npt+i] + (TB[5][1] - 1)*sc[6*npt+i] + (TB[5][2] - 1)*sc[3*npt+i];
        k_f = alpha * k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 7: H + O + M <=> OH + M */
        phi_f = sc[0*npt+i]*sc[1*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[6*npt+i];
        k_f = alpha * k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;

        /*reaction 8: O + O + M <=> O2 + M */
        phi_f = sc[1*npt+i]*sc[1*npt+i];
        alpha = mixture[i] + (TB[7][0] - 1)*sc[14*npt+i] + (TB[7][1] - 1)*sc[4*npt+i] + (TB[7][2] - 1)*sc[6*npt+i];
        k_f = alpha * k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 9: OH + H + M <=> H2O + M */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[8][0] - 1)*sc[3*npt+i] + (TB[8][1] - 1)*sc[6*npt+i];
        k_f = alpha * k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 10: H + O2 <=> O + OH */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[2*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;

        /*reaction 11: H + H + N2 <=> H2 + N2 */
        phi_f = sc[0*npt+i]*sc[0*npt+i]*sc[14*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[14*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 12: H + H + H2 <=> H2 + H2 */
        phi_f = sc[0*npt+i]*sc[0*npt+i]*sc[3*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 13: H + H + H2O <=> H2 + H2O */
        phi_f = sc[0*npt+i]*sc[0*npt+i]*sc[6*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 14: O + H2 <=> OH + H */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[2*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;

        /*reaction 15: O + H2 <=> OH + H */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[2*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;

        /*reaction 16: OH + OH <=> O + H2O */
        phi_f = sc[2*npt+i]*sc[2*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[6*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 17: OH + H2 <=> H + H2O */
        phi_f = sc[2*npt+i]*sc[3*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[6*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 18: H2 + O2 <=> HO2 + H */
        phi_f = sc[3*npt+i]*sc[4*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[5*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 19: HO2 + H <=> OH + OH */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[2*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 20: HO2 + H <=> H2O + O */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[6*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 21: HO2 + O <=> OH + O2 */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 22: HO2 + OH <=> H2O + O2 */
        phi_f = sc[2*npt+i]*sc[5*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 23: HO2 + OH <=> H2O + O2 */
        phi_f = sc[2*npt+i]*sc[5*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 24: HO2 + OH <=> H2O + O2 */
        phi_f = sc[2*npt+i]*sc[5*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[5*npt+i]*sc[5*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[7*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[5*npt+i]*sc[5*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[7*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 27: H2O2 + H <=> H2O + OH */
        phi_f = sc[0*npt+i]*sc[7*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[6*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 28: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[0*npt+i]*sc[7*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 29: H2O2 + O <=> HO2 + OH */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 30: H2O2 + OH <=> H2O + HO2 */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[29*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 31: H2O2 + OH <=> H2O + HO2 */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 32: NO + HO2 <=> NO2 + OH */
        phi_f = sc[5*npt+i]*sc[8*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[9*npt+i];
        Kc = Kc_s[31*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 33: NO2 + H <=> NO + OH */
        phi_f = sc[0*npt+i]*sc[9*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[8*npt+i];
        Kc = Kc_s[32*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;

        /*reaction 34: NO2 + O <=> NO + O2 */
        phi_f = sc[1*npt+i]*sc[9*npt+i];
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[8*npt+i];
        Kc = Kc_s[33*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;

        /*reaction 35: NO2 + NO2 <=> NO + NO + O2 */
        phi_f = sc[9*npt+i]*sc[9*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[8*npt+i]*sc[8*npt+i];
        Kc = Kc_s[34*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;

        /*reaction 36: N2O + H <=> N2 + OH */
        phi_f = sc[0*npt+i]*sc[10*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[14*npt+i];
        Kc = Kc_s[35*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 37: N2O + H <=> N2 + OH */
        phi_f = sc[0*npt+i]*sc[10*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[14*npt+i];
        Kc = Kc_s[36*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 38: N2O + O <=> NO + NO */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[8*npt+i];
        Kc = Kc_s[37*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 39: N2O + O <=> N2 + O2 */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        k_f = k_f_s[38*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[14*npt+i];
        Kc = Kc_s[38*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 40: NH + H <=> N + H2 */
        phi_f = sc[0*npt+i]*sc[11*npt+i];
        k_f = k_f_s[39*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[12*npt+i];
        Kc = Kc_s[39*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 41: NH + O <=> NO + H */
        phi_f = sc[1*npt+i]*sc[11*npt+i];
        k_f = k_f_s[40*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[8*npt+i];
        Kc = Kc_s[40*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 42: NH + OH <=> N + H2O */
        phi_f = sc[2*npt+i]*sc[11*npt+i];
        k_f = k_f_s[41*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[12*npt+i];
        Kc = Kc_s[41*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 43: NH + O2 <=> NO + OH */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[42*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[8*npt+i];
        Kc = Kc_s[42*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 44: NH + NO <=> N2O + H */
        phi_f = sc[8*npt+i]*sc[11*npt+i];
        k_f = k_f_s[43*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[10*npt+i];
        Kc = Kc_s[43*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 45: NH + NO <=> N2O + H */
        phi_f = sc[8*npt+i]*sc[11*npt+i];
        k_f = k_f_s[44*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[10*npt+i];
        Kc = Kc_s[44*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 46: NH + NO <=> N2 + OH */
        phi_f = sc[8*npt+i]*sc[11*npt+i];
        k_f = k_f_s[45*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[14*npt+i];
        Kc = Kc_s[45*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 47: NH + NO2 <=> N2O + OH */
        phi_f = sc[9*npt+i]*sc[11*npt+i];
        k_f = k_f_s[46*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[10*npt+i];
        Kc = Kc_s[46*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 48: N + OH <=> NO + H */
        phi_f = sc[2*npt+i]*sc[12*npt+i];
        k_f = k_f_s[47*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[8*npt+i];
        Kc = Kc_s[47*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 49: N + O2 <=> NO + O */
        phi_f = sc[4*npt+i]*sc[12*npt+i];
        k_f = k_f_s[48*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[8*npt+i];
        Kc = Kc_s[48*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 50: N + NO <=> N2 + O */
        phi_f = sc[8*npt+i]*sc[12*npt+i];
        k_f = k_f_s[49*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[49*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
    }
}

void vcomp_wdot_51_58(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
		double * restrict k_f_s, double * restrict Kc_s,
		double * restrict tc, double * restrict invT, double * restrict T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;

        /*reaction 51: NNH <=> N2 + H */
        phi_f = sc[13*npt+i];
        k_f = k_f_s[50*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[14*npt+i];
        Kc = Kc_s[50*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 52: NNH + H <=> N2 + H2 */
        phi_f = sc[0*npt+i]*sc[13*npt+i];
        k_f = k_f_s[51*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[14*npt+i];
        Kc = Kc_s[51*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 53: NNH + O <=> N2O + H */
        phi_f = sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[52*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[10*npt+i];
        Kc = Kc_s[52*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 54: NNH + O <=> N2 + OH */
        phi_f = sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[53*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[14*npt+i];
        Kc = Kc_s[53*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 55: NNH + O <=> NH + NO */
        phi_f = sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[54*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[11*npt+i];
        Kc = Kc_s[54*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 56: NNH + OH <=> N2 + H2O */
        phi_f = sc[2*npt+i]*sc[13*npt+i];
        k_f = k_f_s[55*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[14*npt+i];
        Kc = Kc_s[55*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 57: NNH + O2 <=> N2 + HO2 */
        phi_f = sc[4*npt+i]*sc[13*npt+i];
        k_f = k_f_s[56*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[14*npt+i];
        Kc = Kc_s[56*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 58: NNH + O2 <=> N2 + H + O2 */
        phi_f = sc[4*npt+i]*sc[13*npt+i];
        k_f = k_f_s[57*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[4*npt+i]*sc[14*npt+i];
        Kc = Kc_s[57*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
    }
}


/*compute the reaction Jacobian */
void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)
{
    double c[15];

    for (int k=0; k<15; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    for (int k=0; k<15; k++) {
        J[240+k] *= 1.e-6;
    }

    /* dTdot/d[X] */
    for (int k=0; k<15; k++) {
        J[k*16+15] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<256; i++) {
        J[i] = 0.0;
    }

    double wdot[15];
    for (int k=0; k<15; k++) {
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
    for (int k = 0; k < 15; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[15];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[15];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[15];
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
    alpha = mixture + (TB[0][0] - 1)*sc[14] + (TB[0][1] - 1)*sc[6] + (TB[0][2] - 1)*sc[3] + (TB[0][3] - 1)*sc[4];
    /* forward */
    phi_f = sc[0]*sc[4];
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
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[0] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[5]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[4] -= q; /* O2 */
    wdot[5] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[4] -= dqdci;                /* dwdot[O2]/d[H] */
        J[5] += dqdci;                /* dwdot[HO2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[0][2] - 1)*dcdc_fac;
        J[48] -= dqdci;               /* dwdot[H]/d[H2] */
        J[52] -= dqdci;               /* dwdot[O2]/d[H2] */
        J[53] += dqdci;               /* dwdot[HO2]/d[H2] */
        /* d()/d[O2] */
        dqdci = (TB[0][3] - 1)*dcdc_fac + k_f*sc[0];
        J[64] -= dqdci;               /* dwdot[H]/d[O2] */
        J[68] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[69] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[80] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[84] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[85] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[96] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[100] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[101] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[N2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[224] -= dqdci;              /* dwdot[H]/d[N2] */
        J[228] -= dqdci;              /* dwdot[O2]/d[N2] */
        J[229] += dqdci;              /* dwdot[HO2]/d[N2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[4];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = TB[0][2]*dcdc_fac;
        dqdc[4] = TB[0][3]*dcdc_fac + k_f*sc[0];
        dqdc[5] = dcdc_fac - k_r;
        dqdc[6] = TB[0][1]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[0][0]*dcdc_fac;
        for (int k=0; k<15; k++) {
            J[16*k+0] -= dqdc[k];
            J[16*k+4] -= dqdc[k];
            J[16*k+5] += dqdc[k];
        }
    }
    J[240] -= dqdT; /* dwdot[H]/dT */
    J[244] -= dqdT; /* dwdot[O2]/dT */
    J[245] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H + O2 (+N2) <=> HO2 (+N2) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = sc[14];
    /* forward */
    phi_f = sc[0]*sc[4];
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
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[0] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[5]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[4] -= q; /* O2 */
    wdot[5] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[4] -= dqdci;                /* dwdot[O2]/d[H] */
        J[5] += dqdci;                /* dwdot[HO2]/d[H] */
        /* d()/d[O2] */
        dqdci =  + k_f*sc[0];
        J[64] -= dqdci;               /* dwdot[H]/d[O2] */
        J[68] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[69] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[80] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[84] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[85] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[N2] */
        dqdci = dcdc_fac;
        J[224] -= dqdci;              /* dwdot[H]/d[N2] */
        J[228] -= dqdci;              /* dwdot[O2]/d[N2] */
        J[229] += dqdci;              /* dwdot[HO2]/d[N2] */
    }
    else {
        dqdc[0] = 0*q_nocor + k_f*sc[4];
        dqdc[1] = 0*q_nocor;
        dqdc[2] = 0*q_nocor;
        dqdc[3] = 0*q_nocor;
        dqdc[4] = 0*q_nocor + k_f*sc[0];
        dqdc[5] = 0*q_nocor - k_r;
        dqdc[6] = 0*q_nocor;
        dqdc[7] = 0*q_nocor;
        dqdc[8] = 0*q_nocor;
        dqdc[9] = 0*q_nocor;
        dqdc[10] = 0*q_nocor;
        dqdc[11] = 0*q_nocor;
        dqdc[12] = 0*q_nocor;
        dqdc[13] = 0*q_nocor;
        dqdc[14] = dcdc_fac;
        for (int k=0; k<15; k++) {
            J[16*k+0] -= dqdc[k];
            J[16*k+4] -= dqdc[k];
            J[16*k+5] += dqdc[k];
        }
    }
    J[240] -= dqdT; /* dwdot[H]/dT */
    J[244] -= dqdT; /* dwdot[O2]/dT */
    J[245] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 3: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[6] + (TB[2][1] - 1)*sc[3];
    /* forward */
    phi_f = sc[7];
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
    phi_r = sc[2]*sc[2];
    Kc = refC * exp(-g_RT[2] - g_RT[2] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2*h_RT[2]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] += 2 * q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[OH] */
        dqdci =  - k_r*2*sc[2];
        J[34] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
        J[39] -= dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2] */
        dqdci = (TB[2][1] - 1)*dcdc_fac;
        J[50] += 2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[55] -= dqdci;               /* dwdot[H2O2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[2][0] - 1)*dcdc_fac;
        J[98] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[103] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
        /* d()/d[H2O2] */
        dqdci =  + k_f;
        J[114] += 2 * dqdci;          /* dwdot[OH]/d[H2O2] */
        J[119] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac - k_r*2*sc[2];
        dqdc[3] = TB[2][1]*dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[2][0]*dcdc_fac;
        dqdc[7] = dcdc_fac + k_f;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        for (int k=0; k<15; k++) {
            J[16*k+2] += 2 * dqdc[k];
            J[16*k+7] -= dqdc[k];
        }
    }
    J[242] += 2 * dqdT; /* dwdot[OH]/dT */
    J[247] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 4: NO + O (+M) <=> NO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[1]*sc[8];
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
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[1] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (h_RT[9]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[8] -= q; /* NO */
    wdot[9] += q; /* NO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[O] */
        dqdci =  + k_f*sc[8];
        J[17] -= dqdci;               /* dwdot[O]/d[O] */
        J[24] -= dqdci;               /* dwdot[NO]/d[O] */
        J[25] += dqdci;               /* dwdot[NO2]/d[O] */
        /* d()/d[NO] */
        dqdci =  + k_f*sc[1];
        J[129] -= dqdci;              /* dwdot[O]/d[NO] */
        J[136] -= dqdci;              /* dwdot[NO]/d[NO] */
        J[137] += dqdci;              /* dwdot[NO2]/d[NO] */
        /* d()/d[NO2] */
        dqdci =  - k_r;
        J[145] -= dqdci;              /* dwdot[O]/d[NO2] */
        J[152] -= dqdci;              /* dwdot[NO]/d[NO2] */
        J[153] += dqdci;              /* dwdot[NO2]/d[NO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[8];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac + k_f*sc[1];
        dqdc[9] = dcdc_fac - k_r;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        for (int k=0; k<15; k++) {
            J[16*k+1] -= dqdc[k];
            J[16*k+8] -= dqdc[k];
            J[16*k+9] += dqdc[k];
        }
    }
    J[241] -= dqdT; /* dwdot[O]/dT */
    J[248] -= dqdT; /* dwdot[NO]/dT */
    J[249] += dqdT; /* dwdot[NO2]/dT */

    /*reaction 5: N2O (+M) <=> N2 + O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[14] + (TB[4][1] - 1)*sc[4] + (TB[4][2] - 1)*sc[6];
    /* forward */
    phi_f = sc[10];
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
    /* Lindemann form */
    F = 1.0;
    dlogFdlogPr = 0.0;
    dlogFdT = 0.0;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = refC * exp(-g_RT[1] + g_RT[10] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10]) + (h_RT[1] + h_RT[14]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[10] -= q; /* N2O */
    wdot[14] += q; /* N2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[O] */
        dqdci =  - k_r*sc[14];
        J[17] += dqdci;               /* dwdot[O]/d[O] */
        J[26] -= dqdci;               /* dwdot[N2O]/d[O] */
        J[30] += dqdci;               /* dwdot[N2]/d[O] */
        /* d()/d[O2] */
        dqdci = (TB[4][1] - 1)*dcdc_fac;
        J[65] += dqdci;               /* dwdot[O]/d[O2] */
        J[74] -= dqdci;               /* dwdot[N2O]/d[O2] */
        J[78] += dqdci;               /* dwdot[N2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[4][2] - 1)*dcdc_fac;
        J[97] += dqdci;               /* dwdot[O]/d[H2O] */
        J[106] -= dqdci;              /* dwdot[N2O]/d[H2O] */
        J[110] += dqdci;              /* dwdot[N2]/d[H2O] */
        /* d()/d[N2O] */
        dqdci =  + k_f;
        J[161] += dqdci;              /* dwdot[O]/d[N2O] */
        J[170] -= dqdci;              /* dwdot[N2O]/d[N2O] */
        J[174] += dqdci;              /* dwdot[N2]/d[N2O] */
        /* d()/d[N2] */
        dqdci = (TB[4][0] - 1)*dcdc_fac - k_r*sc[1];
        J[225] += dqdci;              /* dwdot[O]/d[N2] */
        J[234] -= dqdci;              /* dwdot[N2O]/d[N2] */
        J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac - k_r*sc[14];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[4][1]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[4][2]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac + k_f;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[4][0]*dcdc_fac - k_r*sc[1];
        for (int k=0; k<15; k++) {
            J[16*k+1] += dqdc[k];
            J[16*k+10] -= dqdc[k];
            J[16*k+14] += dqdc[k];
        }
    }
    J[241] += dqdT; /* dwdot[O]/dT */
    J[250] -= dqdT; /* dwdot[N2O]/dT */
    J[254] += dqdT; /* dwdot[N2]/dT */

    /*reaction 6: H + H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[14] + (TB[5][1] - 1)*sc[6] + (TB[5][2] - 1)*sc[3];
    /* forward */
    phi_f = sc[0]*sc[0];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(g_RT[0] + g_RT[0] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[0]) + (h_RT[3]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= 2 * q; /* H */
    wdot[3] += q; /* H2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*2*sc[0];
        J[0] += -2 * dqdci;           /* dwdot[H]/d[H] */
        J[3] += dqdci;                /* dwdot[H2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[5][2] - 1)*q_nocor - k_r;
        J[48] += -2 * dqdci;          /* dwdot[H]/d[H2] */
        J[51] += dqdci;               /* dwdot[H2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*q_nocor;
        J[96] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
        J[99] += dqdci;               /* dwdot[H2]/d[H2O] */
        /* d()/d[N2] */
        dqdci = (TB[5][0] - 1)*q_nocor;
        J[224] += -2 * dqdci;         /* dwdot[H]/d[N2] */
        J[227] += dqdci;              /* dwdot[H2]/d[N2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*2*sc[0];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = TB[5][2] - k_r;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[5][1];
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[5][0];
        for (int k=0; k<15; k++) {
            J[16*k+0] += -2 * dqdc[k];
            J[16*k+3] += dqdc[k];
        }
    }
    J[240] += -2 * dqdT; /* dwdot[H]/dT */
    J[243] += dqdT; /* dwdot[H2]/dT */

    /*reaction 7: H + O + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[6];
    /* forward */
    phi_f = sc[0]*sc[1];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[2];
    Kc = refCinv * exp(g_RT[0] + g_RT[1] - g_RT[2]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[1]) + (h_RT[2]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[1];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[1] -= dqdci;                /* dwdot[O]/d[H] */
        J[2] += dqdci;                /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[0];
        J[16] -= dqdci;               /* dwdot[H]/d[O] */
        J[17] -= dqdci;               /* dwdot[O]/d[O] */
        J[18] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[32] -= dqdci;               /* dwdot[H]/d[OH] */
        J[33] -= dqdci;               /* dwdot[O]/d[OH] */
        J[34] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[6][0] - 1)*q_nocor;
        J[96] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[97] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[98] += dqdci;               /* dwdot[OH]/d[H2O] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[1];
        dqdc[1] = dcdc_fac + k_f*sc[0];
        dqdc[2] = dcdc_fac - k_r;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[6][0];
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        for (int k=0; k<15; k++) {
            J[16*k+0] -= dqdc[k];
            J[16*k+1] -= dqdc[k];
            J[16*k+2] += dqdc[k];
        }
    }
    J[240] -= dqdT; /* dwdot[H]/dT */
    J[241] -= dqdT; /* dwdot[O]/dT */
    J[242] += dqdT; /* dwdot[OH]/dT */

    /*reaction 8: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[7][0] - 1)*sc[14] + (TB[7][1] - 1)*sc[4] + (TB[7][2] - 1)*sc[6];
    /* forward */
    phi_f = sc[1]*sc[1];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[1] + g_RT[1] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[1]) + (h_RT[4]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= 2 * q; /* O */
    wdot[4] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[O] */
        dqdci =  + k_f*2*sc[1];
        J[17] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[20] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[O2] */
        dqdci = (TB[7][1] - 1)*q_nocor - k_r;
        J[65] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        J[68] += dqdci;               /* dwdot[O2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[7][2] - 1)*q_nocor;
        J[97] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        J[100] += dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[N2] */
        dqdci = (TB[7][0] - 1)*q_nocor;
        J[225] += -2 * dqdci;         /* dwdot[O]/d[N2] */
        J[228] += dqdci;              /* dwdot[O2]/d[N2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*2*sc[1];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[7][1] - k_r;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[7][2];
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[7][0];
        for (int k=0; k<15; k++) {
            J[16*k+1] += -2 * dqdc[k];
            J[16*k+4] += dqdc[k];
        }
    }
    J[241] += -2 * dqdT; /* dwdot[O]/dT */
    J[244] += dqdT; /* dwdot[O2]/dT */

    /*reaction 9: OH + H + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[8][0] - 1)*sc[3] + (TB[8][1] - 1)*sc[6];
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] -= q; /* OH */
    wdot[6] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[2] -= dqdci;                /* dwdot[OH]/d[H] */
        J[6] += dqdci;                /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[0];
        J[32] -= dqdci;               /* dwdot[H]/d[OH] */
        J[34] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2] */
        dqdci = (TB[8][0] - 1)*q_nocor;
        J[48] -= dqdci;               /* dwdot[H]/d[H2] */
        J[50] -= dqdci;               /* dwdot[OH]/d[H2] */
        J[54] += dqdci;               /* dwdot[H2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[8][1] - 1)*q_nocor - k_r;
        J[96] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[98] -= dqdci;               /* dwdot[OH]/d[H2O] */
        J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[2];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[0];
        dqdc[3] = TB[8][0];
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[8][1] - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        for (int k=0; k<15; k++) {
            J[16*k+0] -= dqdc[k];
            J[16*k+2] -= dqdc[k];
            J[16*k+6] += dqdc[k];
        }
    }
    J[240] -= dqdT; /* dwdot[H]/dT */
    J[242] -= dqdT; /* dwdot[OH]/dT */
    J[246] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 10: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[2];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[2] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[1] + h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* O */
    wdot[2] += q; /* OH */
    wdot[4] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[O]/d[H] */
    J[2] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] -= dqdci;                /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[2];
    J[16] -= dqdci;               /* dwdot[H]/d[O] */
    J[17] += dqdci;               /* dwdot[O]/d[O] */
    J[18] += dqdci;               /* dwdot[OH]/d[O] */
    J[20] -= dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[32] -= dqdci;               /* dwdot[H]/d[OH] */
    J[33] += dqdci;               /* dwdot[O]/d[OH] */
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[36] -= dqdci;               /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[64] -= dqdci;               /* dwdot[H]/d[O2] */
    J[65] += dqdci;               /* dwdot[O]/d[O2] */
    J[66] += dqdci;               /* dwdot[OH]/d[O2] */
    J[68] -= dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[240] -= dqdT;               /* dwdot[H]/dT */
    J[241] += dqdT;               /* dwdot[O]/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[244] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 11: H + H + N2 <=> H2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[0]*sc[14];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[14];
    Kc = refCinv * exp(g_RT[0] + g_RT[0] - g_RT[3] + g_RT[14] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[0] + h_RT[14]) + (h_RT[3] + h_RT[14]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= 2 * q; /* H */
    wdot[3] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*2*sc[0]*sc[14];
    J[0] += -2 * dqdci;           /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[14];
    J[48] += -2 * dqdci;          /* dwdot[H]/d[H2] */
    J[51] += dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[0]*sc[0] - k_r*sc[3];
    J[224] += -2 * dqdci;         /* dwdot[H]/d[N2] */
    J[227] += dqdci;              /* dwdot[H2]/d[N2] */
    /* d()/dT */
    J[240] += -2 * dqdT;          /* dwdot[H]/dT */
    J[243] += dqdT;               /* dwdot[H2]/dT */

    /*reaction 12: H + H + H2 <=> H2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[0]*sc[3];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[3];
    Kc = refCinv * exp(g_RT[0] + g_RT[0] + g_RT[3] - g_RT[3] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[0] + h_RT[3]) + (2*h_RT[3]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= 2 * q; /* H */
    wdot[3] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*2*sc[0]*sc[3];
    J[0] += -2 * dqdci;           /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[0]*sc[0] - k_r*2*sc[3];
    J[48] += -2 * dqdci;          /* dwdot[H]/d[H2] */
    J[51] += dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[240] += -2 * dqdT;          /* dwdot[H]/dT */
    J[243] += dqdT;               /* dwdot[H2]/dT */

    /*reaction 13: H + H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[0]*sc[6];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = refCinv * exp(g_RT[0] + g_RT[0] - g_RT[3] + g_RT[6] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[0] + h_RT[6]) + (h_RT[3] + h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= 2 * q; /* H */
    wdot[3] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*2*sc[0]*sc[6];
    J[0] += -2 * dqdci;           /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[48] += -2 * dqdci;          /* dwdot[H]/d[H2] */
    J[51] += dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[0]*sc[0] - k_r*sc[3];
    J[96] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
    J[99] += dqdci;               /* dwdot[H2]/d[H2O] */
    /* d()/dT */
    J[240] += -2 * dqdT;          /* dwdot[H]/dT */
    J[243] += dqdT;               /* dwdot[H2]/dT */

    /*reaction 14: O + H2 <=> OH + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[2];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[2] + g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* OH */
    wdot[3] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[2];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[O]/d[H] */
    J[2] += dqdci;                /* dwdot[OH]/d[H] */
    J[3] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[3];
    J[16] += dqdci;               /* dwdot[H]/d[O] */
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[18] += dqdci;               /* dwdot[OH]/d[O] */
    J[19] -= dqdci;               /* dwdot[H2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[32] += dqdci;               /* dwdot[H]/d[OH] */
    J[33] -= dqdci;               /* dwdot[O]/d[OH] */
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[35] -= dqdci;               /* dwdot[H2]/d[OH] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[1];
    J[48] += dqdci;               /* dwdot[H]/d[H2] */
    J[49] -= dqdci;               /* dwdot[O]/d[H2] */
    J[50] += dqdci;               /* dwdot[OH]/d[H2] */
    J[51] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[243] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 15: O + H2 <=> OH + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[2];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[2] + g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* OH */
    wdot[3] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[2];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[O]/d[H] */
    J[2] += dqdci;                /* dwdot[OH]/d[H] */
    J[3] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[3];
    J[16] += dqdci;               /* dwdot[H]/d[O] */
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[18] += dqdci;               /* dwdot[OH]/d[O] */
    J[19] -= dqdci;               /* dwdot[H2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[32] += dqdci;               /* dwdot[H]/d[OH] */
    J[33] -= dqdci;               /* dwdot[O]/d[OH] */
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[35] -= dqdci;               /* dwdot[H2]/d[OH] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[1];
    J[48] += dqdci;               /* dwdot[H]/d[H2] */
    J[49] -= dqdci;               /* dwdot[O]/d[H2] */
    J[50] += dqdci;               /* dwdot[OH]/d[H2] */
    J[51] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[243] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 16: OH + OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[2];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[6];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[2] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[2]) + (h_RT[1] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[2] -= 2 * q; /* OH */
    wdot[6] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[6];
    J[17] += dqdci;               /* dwdot[O]/d[O] */
    J[18] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[22] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2*sc[2];
    J[33] += dqdci;               /* dwdot[O]/d[OH] */
    J[34] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[97] += dqdci;               /* dwdot[O]/d[H2O] */
    J[98] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[241] += dqdT;               /* dwdot[O]/dT */
    J[242] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 17: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[2] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[2] -= q; /* OH */
    wdot[3] -= q; /* H2 */
    wdot[6] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[2] -= dqdci;                /* dwdot[OH]/d[H] */
    J[3] -= dqdci;                /* dwdot[H2]/d[H] */
    J[6] += dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[3];
    J[32] += dqdci;               /* dwdot[H]/d[OH] */
    J[34] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[35] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[48] += dqdci;               /* dwdot[H]/d[H2] */
    J[50] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[51] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[54] += dqdci;               /* dwdot[H2O]/d[H2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0];
    J[96] += dqdci;               /* dwdot[H]/d[H2O] */
    J[98] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[99] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[242] -= dqdT;               /* dwdot[OH]/dT */
    J[243] -= dqdT;               /* dwdot[H2]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 18: H2 + O2 <=> HO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[4];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[5];
    Kc = exp(-g_RT[0] + g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[4]) + (h_RT[0] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] -= q; /* H2 */
    wdot[4] -= q; /* O2 */
    wdot[5] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] -= dqdci;                /* dwdot[H2]/d[H] */
    J[4] -= dqdci;                /* dwdot[O2]/d[H] */
    J[5] += dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[48] += dqdci;               /* dwdot[H]/d[H2] */
    J[51] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[52] -= dqdci;               /* dwdot[O2]/d[H2] */
    J[53] += dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[3];
    J[64] += dqdci;               /* dwdot[H]/d[O2] */
    J[67] -= dqdci;               /* dwdot[H2]/d[O2] */
    J[68] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[69] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[80] += dqdci;               /* dwdot[H]/d[HO2] */
    J[83] -= dqdci;               /* dwdot[H2]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[85] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[243] -= dqdT;               /* dwdot[H2]/dT */
    J[244] -= dqdT;               /* dwdot[O2]/dT */
    J[245] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[2];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[2] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (2*h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += 2 * q; /* OH */
    wdot[5] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += 2 * dqdci;            /* dwdot[OH]/d[H] */
    J[5] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[2];
    J[32] -= dqdci;               /* dwdot[H]/d[OH] */
    J[34] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[37] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[80] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[82] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[240] -= dqdT;               /* dwdot[H]/dT */
    J[242] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[245] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 20: HO2 + H <=> H2O + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
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
    wdot[0] -= q; /* H */
    wdot[1] += q; /* O */
    wdot[5] -= q; /* HO2 */
    wdot[6] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[O]/d[H] */
    J[5] -= dqdci;                /* dwdot[HO2]/d[H] */
    J[6] += dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[6];
    J[16] -= dqdci;               /* dwdot[H]/d[O] */
    J[17] += dqdci;               /* dwdot[O]/d[O] */
    J[21] -= dqdci;               /* dwdot[HO2]/d[O] */
    J[22] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[80] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[81] += dqdci;               /* dwdot[O]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    J[86] += dqdci;               /* dwdot[H2O]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[96] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[97] += dqdci;               /* dwdot[O]/d[H2O] */
    J[101] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[240] -= dqdT;               /* dwdot[H]/dT */
    J[241] += dqdT;               /* dwdot[O]/dT */
    J[245] -= dqdT;               /* dwdot[HO2]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 21: HO2 + O <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[4] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* OH */
    wdot[4] += q; /* O2 */
    wdot[5] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[5];
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[18] += dqdci;               /* dwdot[OH]/d[O] */
    J[20] += dqdci;               /* dwdot[O2]/d[O] */
    J[21] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[33] -= dqdci;               /* dwdot[O]/d[OH] */
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[36] += dqdci;               /* dwdot[O2]/d[OH] */
    J[37] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[2];
    J[65] -= dqdci;               /* dwdot[O]/d[O2] */
    J[66] += dqdci;               /* dwdot[OH]/d[O2] */
    J[68] += dqdci;               /* dwdot[O2]/d[O2] */
    J[69] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[81] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[82] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[84] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[244] += dqdT;               /* dwdot[O2]/dT */
    J[245] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 22: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[2] - g_RT[4] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* OH */
    wdot[4] += q; /* O2 */
    wdot[5] -= q; /* HO2 */
    wdot[6] += q; /* H2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[5];
    J[34] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[36] += dqdci;               /* dwdot[O2]/d[OH] */
    J[37] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[6];
    J[66] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[68] += dqdci;               /* dwdot[O2]/d[O2] */
    J[69] -= dqdci;               /* dwdot[HO2]/d[O2] */
    J[70] += dqdci;               /* dwdot[H2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[82] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[84] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    J[86] += dqdci;               /* dwdot[H2O]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[98] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[100] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[101] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[242] -= dqdT;               /* dwdot[OH]/dT */
    J[244] += dqdT;               /* dwdot[O2]/dT */
    J[245] -= dqdT;               /* dwdot[HO2]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[2] - g_RT[4] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* OH */
    wdot[4] += q; /* O2 */
    wdot[5] -= q; /* HO2 */
    wdot[6] += q; /* H2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[5];
    J[34] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[36] += dqdci;               /* dwdot[O2]/d[OH] */
    J[37] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[6];
    J[66] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[68] += dqdci;               /* dwdot[O2]/d[O2] */
    J[69] -= dqdci;               /* dwdot[HO2]/d[O2] */
    J[70] += dqdci;               /* dwdot[H2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[82] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[84] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    J[86] += dqdci;               /* dwdot[H2O]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[98] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[100] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[101] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[242] -= dqdT;               /* dwdot[OH]/dT */
    J[244] += dqdT;               /* dwdot[O2]/dT */
    J[245] -= dqdT;               /* dwdot[HO2]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[2] - g_RT[4] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* OH */
    wdot[4] += q; /* O2 */
    wdot[5] -= q; /* HO2 */
    wdot[6] += q; /* H2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[5];
    J[34] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[36] += dqdci;               /* dwdot[O2]/d[OH] */
    J[37] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[6];
    J[66] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[68] += dqdci;               /* dwdot[O2]/d[O2] */
    J[69] -= dqdci;               /* dwdot[HO2]/d[O2] */
    J[70] += dqdci;               /* dwdot[H2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[82] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[84] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    J[86] += dqdci;               /* dwdot[H2O]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[98] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[100] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[101] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[242] -= dqdT;               /* dwdot[OH]/dT */
    J[244] += dqdT;               /* dwdot[O2]/dT */
    J[245] -= dqdT;               /* dwdot[HO2]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[5];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[7];
    Kc = exp(-g_RT[4] + g_RT[5] + g_RT[5] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[5]) + (h_RT[4] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* O2 */
    wdot[5] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[68] += dqdci;               /* dwdot[O2]/d[O2] */
    J[69] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[5];
    J[84] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[85] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[87] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[4];
    J[116] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[117] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[119] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[244] += dqdT;               /* dwdot[O2]/dT */
    J[245] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[247] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[5];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[7];
    Kc = exp(-g_RT[4] + g_RT[5] + g_RT[5] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[5]) + (h_RT[4] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* O2 */
    wdot[5] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[68] += dqdci;               /* dwdot[O2]/d[O2] */
    J[69] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[5];
    J[84] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[85] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[87] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[4];
    J[116] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[117] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[119] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[244] += dqdT;               /* dwdot[O2]/dT */
    J[245] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[247] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 27: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[6];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[OH]/d[H] */
    J[6] += dqdci;                /* dwdot[H2O]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[32] -= dqdci;               /* dwdot[H]/d[OH] */
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[39] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[96] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[98] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[103] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[112] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[114] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[118] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[119] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[240] -= dqdT;               /* dwdot[H]/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */
    J[247] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 28: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* H2 */
    wdot[5] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[H2]/d[H] */
    J[5] += dqdci;                /* dwdot[HO2]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[48] -= dqdci;               /* dwdot[H]/d[H2] */
    J[51] += dqdci;               /* dwdot[H2]/d[H2] */
    J[53] += dqdci;               /* dwdot[HO2]/d[H2] */
    J[55] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[80] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[83] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[85] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[87] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[112] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[115] += dqdci;              /* dwdot[H2]/d[H2O2] */
    J[117] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[119] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[240] -= dqdT;               /* dwdot[H]/dT */
    J[243] += dqdT;               /* dwdot[H2]/dT */
    J[245] += dqdT;               /* dwdot[HO2]/dT */
    J[247] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* OH */
    wdot[5] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[18] += dqdci;               /* dwdot[OH]/d[O] */
    J[21] += dqdci;               /* dwdot[HO2]/d[O] */
    J[23] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[33] -= dqdci;               /* dwdot[O]/d[OH] */
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[37] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[39] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[2];
    J[81] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[82] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[85] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[87] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[113] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[114] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[117] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[119] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[245] += dqdT;               /* dwdot[HO2]/dT */
    J[247] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 30: H2O2 + OH <=> H2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[2] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* OH */
    wdot[5] += q; /* HO2 */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[34] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[37] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[39] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[6];
    J[82] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[85] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[86] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[87] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[98] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[101] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[103] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[117] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[118] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[119] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[242] -= dqdT;               /* dwdot[OH]/dT */
    J[245] += dqdT;               /* dwdot[HO2]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */
    J[247] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 31: H2O2 + OH <=> H2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[2] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* OH */
    wdot[5] += q; /* HO2 */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[34] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[37] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[39] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[6];
    J[82] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[85] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[86] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[87] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[98] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[101] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[103] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[117] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[118] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[119] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[242] -= dqdT;               /* dwdot[OH]/dT */
    J[245] += dqdT;               /* dwdot[HO2]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */
    J[247] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 32: NO + HO2 <=> NO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[8];
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9];
    Kc = exp(-g_RT[2] + g_RT[5] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[8]) + (h_RT[2] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* OH */
    wdot[5] -= q; /* HO2 */
    wdot[8] -= q; /* NO */
    wdot[9] += q; /* NO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[37] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[40] -= dqdci;               /* dwdot[NO]/d[OH] */
    J[41] += dqdci;               /* dwdot[NO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[8];
    J[82] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    J[88] -= dqdci;               /* dwdot[NO]/d[HO2] */
    J[89] += dqdci;               /* dwdot[NO2]/d[HO2] */
    /* d()/d[NO] */
    dqdci =  + k_f*sc[5];
    J[130] += dqdci;              /* dwdot[OH]/d[NO] */
    J[133] -= dqdci;              /* dwdot[HO2]/d[NO] */
    J[136] -= dqdci;              /* dwdot[NO]/d[NO] */
    J[137] += dqdci;              /* dwdot[NO2]/d[NO] */
    /* d()/d[NO2] */
    dqdci =  - k_r*sc[2];
    J[146] += dqdci;              /* dwdot[OH]/d[NO2] */
    J[149] -= dqdci;              /* dwdot[HO2]/d[NO2] */
    J[152] -= dqdci;              /* dwdot[NO]/d[NO2] */
    J[153] += dqdci;              /* dwdot[NO2]/d[NO2] */
    /* d()/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[245] -= dqdT;               /* dwdot[HO2]/dT */
    J[248] -= dqdT;               /* dwdot[NO]/dT */
    J[249] += dqdT;               /* dwdot[NO2]/dT */

    /*reaction 33: NO2 + H <=> NO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[9];
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[8];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[8] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[9]) + (h_RT[2] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* OH */
    wdot[8] += q; /* NO */
    wdot[9] -= q; /* NO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[9];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[OH]/d[H] */
    J[8] += dqdci;                /* dwdot[NO]/d[H] */
    J[9] -= dqdci;                /* dwdot[NO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[8];
    J[32] -= dqdci;               /* dwdot[H]/d[OH] */
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[40] += dqdci;               /* dwdot[NO]/d[OH] */
    J[41] -= dqdci;               /* dwdot[NO2]/d[OH] */
    /* d()/d[NO] */
    dqdci =  - k_r*sc[2];
    J[128] -= dqdci;              /* dwdot[H]/d[NO] */
    J[130] += dqdci;              /* dwdot[OH]/d[NO] */
    J[136] += dqdci;              /* dwdot[NO]/d[NO] */
    J[137] -= dqdci;              /* dwdot[NO2]/d[NO] */
    /* d()/d[NO2] */
    dqdci =  + k_f*sc[0];
    J[144] -= dqdci;              /* dwdot[H]/d[NO2] */
    J[146] += dqdci;              /* dwdot[OH]/d[NO2] */
    J[152] += dqdci;              /* dwdot[NO]/d[NO2] */
    J[153] -= dqdci;              /* dwdot[NO2]/d[NO2] */
    /* d()/dT */
    J[240] -= dqdT;               /* dwdot[H]/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[248] += dqdT;               /* dwdot[NO]/dT */
    J[249] -= dqdT;               /* dwdot[NO2]/dT */

    /*reaction 34: NO2 + O <=> NO + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[8];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[8] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (h_RT[4] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[4] += q; /* O2 */
    wdot[8] += q; /* NO */
    wdot[9] -= q; /* NO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[20] += dqdci;               /* dwdot[O2]/d[O] */
    J[24] += dqdci;               /* dwdot[NO]/d[O] */
    J[25] -= dqdci;               /* dwdot[NO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[8];
    J[65] -= dqdci;               /* dwdot[O]/d[O2] */
    J[68] += dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[NO]/d[O2] */
    J[73] -= dqdci;               /* dwdot[NO2]/d[O2] */
    /* d()/d[NO] */
    dqdci =  - k_r*sc[4];
    J[129] -= dqdci;              /* dwdot[O]/d[NO] */
    J[132] += dqdci;              /* dwdot[O2]/d[NO] */
    J[136] += dqdci;              /* dwdot[NO]/d[NO] */
    J[137] -= dqdci;              /* dwdot[NO2]/d[NO] */
    /* d()/d[NO2] */
    dqdci =  + k_f*sc[1];
    J[145] -= dqdci;              /* dwdot[O]/d[NO2] */
    J[148] += dqdci;              /* dwdot[O2]/d[NO2] */
    J[152] += dqdci;              /* dwdot[NO]/d[NO2] */
    J[153] -= dqdci;              /* dwdot[NO2]/d[NO2] */
    /* d()/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[244] += dqdT;               /* dwdot[O2]/dT */
    J[248] += dqdT;               /* dwdot[NO]/dT */
    J[249] -= dqdT;               /* dwdot[NO2]/dT */

    /*reaction 35: NO2 + NO2 <=> NO + NO + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[9];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[8]*sc[8];
    Kc = refC * exp(-g_RT[4] - g_RT[8] - g_RT[8] + g_RT[9] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[9]) + (h_RT[4] + 2*h_RT[8]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* O2 */
    wdot[8] += 2 * q; /* NO */
    wdot[9] -= 2 * q; /* NO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[8]*sc[8];
    J[68] += dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += 2 * dqdci;           /* dwdot[NO]/d[O2] */
    J[73] += -2 * dqdci;          /* dwdot[NO2]/d[O2] */
    /* d()/d[NO] */
    dqdci =  - k_r*sc[4]*2*sc[8];
    J[132] += dqdci;              /* dwdot[O2]/d[NO] */
    J[136] += 2 * dqdci;          /* dwdot[NO]/d[NO] */
    J[137] += -2 * dqdci;         /* dwdot[NO2]/d[NO] */
    /* d()/d[NO2] */
    dqdci =  + k_f*2*sc[9];
    J[148] += dqdci;              /* dwdot[O2]/d[NO2] */
    J[152] += 2 * dqdci;          /* dwdot[NO]/d[NO2] */
    J[153] += -2 * dqdci;         /* dwdot[NO2]/d[NO2] */
    /* d()/dT */
    J[244] += dqdT;               /* dwdot[O2]/dT */
    J[248] += 2 * dqdT;           /* dwdot[NO]/dT */
    J[249] += -2 * dqdT;          /* dwdot[NO2]/dT */

    /*reaction 36: N2O + H <=> N2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[10];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[14];
    Kc = exp(g_RT[0] - g_RT[2] + g_RT[10] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[10]) + (h_RT[2] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* OH */
    wdot[10] -= q; /* N2O */
    wdot[14] += q; /* N2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[OH]/d[H] */
    J[10] -= dqdci;               /* dwdot[N2O]/d[H] */
    J[14] += dqdci;               /* dwdot[N2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[32] -= dqdci;               /* dwdot[H]/d[OH] */
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[42] -= dqdci;               /* dwdot[N2O]/d[OH] */
    J[46] += dqdci;               /* dwdot[N2]/d[OH] */
    /* d()/d[N2O] */
    dqdci =  + k_f*sc[0];
    J[160] -= dqdci;              /* dwdot[H]/d[N2O] */
    J[162] += dqdci;              /* dwdot[OH]/d[N2O] */
    J[170] -= dqdci;              /* dwdot[N2O]/d[N2O] */
    J[174] += dqdci;              /* dwdot[N2]/d[N2O] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[2];
    J[224] -= dqdci;              /* dwdot[H]/d[N2] */
    J[226] += dqdci;              /* dwdot[OH]/d[N2] */
    J[234] -= dqdci;              /* dwdot[N2O]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[240] -= dqdT;               /* dwdot[H]/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[250] -= dqdT;               /* dwdot[N2O]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[10];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[14];
    Kc = exp(g_RT[0] - g_RT[2] + g_RT[10] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[10]) + (h_RT[2] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* OH */
    wdot[10] -= q; /* N2O */
    wdot[14] += q; /* N2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[OH]/d[H] */
    J[10] -= dqdci;               /* dwdot[N2O]/d[H] */
    J[14] += dqdci;               /* dwdot[N2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[32] -= dqdci;               /* dwdot[H]/d[OH] */
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[42] -= dqdci;               /* dwdot[N2O]/d[OH] */
    J[46] += dqdci;               /* dwdot[N2]/d[OH] */
    /* d()/d[N2O] */
    dqdci =  + k_f*sc[0];
    J[160] -= dqdci;              /* dwdot[H]/d[N2O] */
    J[162] += dqdci;              /* dwdot[OH]/d[N2O] */
    J[170] -= dqdci;              /* dwdot[N2O]/d[N2O] */
    J[174] += dqdci;              /* dwdot[N2]/d[N2O] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[2];
    J[224] -= dqdci;              /* dwdot[H]/d[N2] */
    J[226] += dqdci;              /* dwdot[OH]/d[N2] */
    J[234] -= dqdci;              /* dwdot[N2O]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[240] -= dqdT;               /* dwdot[H]/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[250] -= dqdT;               /* dwdot[N2O]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 38: N2O + O <=> NO + NO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[8];
    Kc = exp(g_RT[1] - g_RT[8] - g_RT[8] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (2*h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[8] += 2 * q; /* NO */
    wdot[10] -= q; /* N2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += 2 * dqdci;           /* dwdot[NO]/d[O] */
    J[26] -= dqdci;               /* dwdot[N2O]/d[O] */
    /* d()/d[NO] */
    dqdci =  - k_r*2*sc[8];
    J[129] -= dqdci;              /* dwdot[O]/d[NO] */
    J[136] += 2 * dqdci;          /* dwdot[NO]/d[NO] */
    J[138] -= dqdci;              /* dwdot[N2O]/d[NO] */
    /* d()/d[N2O] */
    dqdci =  + k_f*sc[1];
    J[161] -= dqdci;              /* dwdot[O]/d[N2O] */
    J[168] += 2 * dqdci;          /* dwdot[NO]/d[N2O] */
    J[170] -= dqdci;              /* dwdot[N2O]/d[N2O] */
    /* d()/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[248] += 2 * dqdT;           /* dwdot[NO]/dT */
    J[250] -= dqdT;               /* dwdot[N2O]/dT */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = prefactor_units[38] * fwd_A[38]
                * exp(fwd_beta[38] * tc[0] - activation_units[38] * fwd_Ea[38] * invT);
    dlnkfdT = fwd_beta[38] * invT + activation_units[38] * fwd_Ea[38] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(g_RT[1] - g_RT[4] + g_RT[10] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[4] += q; /* O2 */
    wdot[10] -= q; /* N2O */
    wdot[14] += q; /* N2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[20] += dqdci;               /* dwdot[O2]/d[O] */
    J[26] -= dqdci;               /* dwdot[N2O]/d[O] */
    J[30] += dqdci;               /* dwdot[N2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[14];
    J[65] -= dqdci;               /* dwdot[O]/d[O2] */
    J[68] += dqdci;               /* dwdot[O2]/d[O2] */
    J[74] -= dqdci;               /* dwdot[N2O]/d[O2] */
    J[78] += dqdci;               /* dwdot[N2]/d[O2] */
    /* d()/d[N2O] */
    dqdci =  + k_f*sc[1];
    J[161] -= dqdci;              /* dwdot[O]/d[N2O] */
    J[164] += dqdci;              /* dwdot[O2]/d[N2O] */
    J[170] -= dqdci;              /* dwdot[N2O]/d[N2O] */
    J[174] += dqdci;              /* dwdot[N2]/d[N2O] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[4];
    J[225] -= dqdci;              /* dwdot[O]/d[N2] */
    J[228] += dqdci;              /* dwdot[O2]/d[N2] */
    J[234] -= dqdci;              /* dwdot[N2O]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[244] += dqdT;               /* dwdot[O2]/dT */
    J[250] -= dqdT;               /* dwdot[N2O]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 40: NH + H <=> N + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[11];
    k_f = prefactor_units[39] * fwd_A[39]
                * exp(fwd_beta[39] * tc[0] - activation_units[39] * fwd_Ea[39] * invT);
    dlnkfdT = fwd_beta[39] * invT + activation_units[39] * fwd_Ea[39] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[12];
    Kc = exp(g_RT[0] - g_RT[3] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[11]) + (h_RT[3] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* H2 */
    wdot[11] -= q; /* NH */
    wdot[12] += q; /* N */
    /* d()/d[H] */
    dqdci =  + k_f*sc[11];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[H2]/d[H] */
    J[11] -= dqdci;               /* dwdot[NH]/d[H] */
    J[12] += dqdci;               /* dwdot[N]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[12];
    J[48] -= dqdci;               /* dwdot[H]/d[H2] */
    J[51] += dqdci;               /* dwdot[H2]/d[H2] */
    J[59] -= dqdci;               /* dwdot[NH]/d[H2] */
    J[60] += dqdci;               /* dwdot[N]/d[H2] */
    /* d()/d[NH] */
    dqdci =  + k_f*sc[0];
    J[176] -= dqdci;              /* dwdot[H]/d[NH] */
    J[179] += dqdci;              /* dwdot[H2]/d[NH] */
    J[187] -= dqdci;              /* dwdot[NH]/d[NH] */
    J[188] += dqdci;              /* dwdot[N]/d[NH] */
    /* d()/d[N] */
    dqdci =  - k_r*sc[3];
    J[192] -= dqdci;              /* dwdot[H]/d[N] */
    J[195] += dqdci;              /* dwdot[H2]/d[N] */
    J[203] -= dqdci;              /* dwdot[NH]/d[N] */
    J[204] += dqdci;              /* dwdot[N]/d[N] */
    /* d()/dT */
    J[240] -= dqdT;               /* dwdot[H]/dT */
    J[243] += dqdT;               /* dwdot[H2]/dT */
    J[251] -= dqdT;               /* dwdot[NH]/dT */
    J[252] += dqdT;               /* dwdot[N]/dT */

    /*reaction 41: NH + O <=> NO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = prefactor_units[40] * fwd_A[40]
                * exp(fwd_beta[40] * tc[0] - activation_units[40] * fwd_Ea[40] * invT);
    dlnkfdT = fwd_beta[40] * invT + activation_units[40] * fwd_Ea[40] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[8];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[8] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[11]) + (h_RT[0] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* O */
    wdot[8] += q; /* NO */
    wdot[11] -= q; /* NH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[8];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[O]/d[H] */
    J[8] += dqdci;                /* dwdot[NO]/d[H] */
    J[11] -= dqdci;               /* dwdot[NH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[11];
    J[16] += dqdci;               /* dwdot[H]/d[O] */
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += dqdci;               /* dwdot[NO]/d[O] */
    J[27] -= dqdci;               /* dwdot[NH]/d[O] */
    /* d()/d[NO] */
    dqdci =  - k_r*sc[0];
    J[128] += dqdci;              /* dwdot[H]/d[NO] */
    J[129] -= dqdci;              /* dwdot[O]/d[NO] */
    J[136] += dqdci;              /* dwdot[NO]/d[NO] */
    J[139] -= dqdci;              /* dwdot[NH]/d[NO] */
    /* d()/d[NH] */
    dqdci =  + k_f*sc[1];
    J[176] += dqdci;              /* dwdot[H]/d[NH] */
    J[177] -= dqdci;              /* dwdot[O]/d[NH] */
    J[184] += dqdci;              /* dwdot[NO]/d[NH] */
    J[187] -= dqdci;              /* dwdot[NH]/d[NH] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[248] += dqdT;               /* dwdot[NO]/dT */
    J[251] -= dqdT;               /* dwdot[NH]/dT */

    /*reaction 42: NH + OH <=> N + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = prefactor_units[41] * fwd_A[41]
                * exp(fwd_beta[41] * tc[0] - activation_units[41] * fwd_Ea[41] * invT);
    dlnkfdT = fwd_beta[41] * invT + activation_units[41] * fwd_Ea[41] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[12];
    Kc = exp(g_RT[2] - g_RT[6] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[6] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[11] -= q; /* NH */
    wdot[12] += q; /* N */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[34] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[43] -= dqdci;               /* dwdot[NH]/d[OH] */
    J[44] += dqdci;               /* dwdot[N]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[12];
    J[98] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[107] -= dqdci;              /* dwdot[NH]/d[H2O] */
    J[108] += dqdci;              /* dwdot[N]/d[H2O] */
    /* d()/d[NH] */
    dqdci =  + k_f*sc[2];
    J[178] -= dqdci;              /* dwdot[OH]/d[NH] */
    J[182] += dqdci;              /* dwdot[H2O]/d[NH] */
    J[187] -= dqdci;              /* dwdot[NH]/d[NH] */
    J[188] += dqdci;              /* dwdot[N]/d[NH] */
    /* d()/d[N] */
    dqdci =  - k_r*sc[6];
    J[194] -= dqdci;              /* dwdot[OH]/d[N] */
    J[198] += dqdci;              /* dwdot[H2O]/d[N] */
    J[203] -= dqdci;              /* dwdot[NH]/d[N] */
    J[204] += dqdci;              /* dwdot[N]/d[N] */
    /* d()/dT */
    J[242] -= dqdT;               /* dwdot[OH]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */
    J[251] -= dqdT;               /* dwdot[NH]/dT */
    J[252] += dqdT;               /* dwdot[N]/dT */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[42] * fwd_A[42]
                * exp(fwd_beta[42] * tc[0] - activation_units[42] * fwd_Ea[42] * invT);
    dlnkfdT = fwd_beta[42] * invT + activation_units[42] * fwd_Ea[42] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[8];
    Kc = exp(-g_RT[2] + g_RT[4] - g_RT[8] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[2] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* OH */
    wdot[4] -= q; /* O2 */
    wdot[8] += q; /* NO */
    wdot[11] -= q; /* NH */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[8];
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[36] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[40] += dqdci;               /* dwdot[NO]/d[OH] */
    J[43] -= dqdci;               /* dwdot[NH]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[66] += dqdci;               /* dwdot[OH]/d[O2] */
    J[68] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[NO]/d[O2] */
    J[75] -= dqdci;               /* dwdot[NH]/d[O2] */
    /* d()/d[NO] */
    dqdci =  - k_r*sc[2];
    J[130] += dqdci;              /* dwdot[OH]/d[NO] */
    J[132] -= dqdci;              /* dwdot[O2]/d[NO] */
    J[136] += dqdci;              /* dwdot[NO]/d[NO] */
    J[139] -= dqdci;              /* dwdot[NH]/d[NO] */
    /* d()/d[NH] */
    dqdci =  + k_f*sc[4];
    J[178] += dqdci;              /* dwdot[OH]/d[NH] */
    J[180] -= dqdci;              /* dwdot[O2]/d[NH] */
    J[184] += dqdci;              /* dwdot[NO]/d[NH] */
    J[187] -= dqdci;              /* dwdot[NH]/d[NH] */
    /* d()/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[244] -= dqdT;               /* dwdot[O2]/dT */
    J[248] += dqdT;               /* dwdot[NO]/dT */
    J[251] -= dqdT;               /* dwdot[NH]/dT */

    /*reaction 44: NH + NO <=> N2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = prefactor_units[43] * fwd_A[43]
                * exp(fwd_beta[43] * tc[0] - activation_units[43] * fwd_Ea[43] * invT);
    dlnkfdT = fwd_beta[43] * invT + activation_units[43] * fwd_Ea[43] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[10];
    Kc = exp(-g_RT[0] + g_RT[8] - g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (h_RT[0] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[8] -= q; /* NO */
    wdot[10] += q; /* N2O */
    wdot[11] -= q; /* NH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[10];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[8] -= dqdci;                /* dwdot[NO]/d[H] */
    J[10] += dqdci;               /* dwdot[N2O]/d[H] */
    J[11] -= dqdci;               /* dwdot[NH]/d[H] */
    /* d()/d[NO] */
    dqdci =  + k_f*sc[11];
    J[128] += dqdci;              /* dwdot[H]/d[NO] */
    J[136] -= dqdci;              /* dwdot[NO]/d[NO] */
    J[138] += dqdci;              /* dwdot[N2O]/d[NO] */
    J[139] -= dqdci;              /* dwdot[NH]/d[NO] */
    /* d()/d[N2O] */
    dqdci =  - k_r*sc[0];
    J[160] += dqdci;              /* dwdot[H]/d[N2O] */
    J[168] -= dqdci;              /* dwdot[NO]/d[N2O] */
    J[170] += dqdci;              /* dwdot[N2O]/d[N2O] */
    J[171] -= dqdci;              /* dwdot[NH]/d[N2O] */
    /* d()/d[NH] */
    dqdci =  + k_f*sc[8];
    J[176] += dqdci;              /* dwdot[H]/d[NH] */
    J[184] -= dqdci;              /* dwdot[NO]/d[NH] */
    J[186] += dqdci;              /* dwdot[N2O]/d[NH] */
    J[187] -= dqdci;              /* dwdot[NH]/d[NH] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[248] -= dqdT;               /* dwdot[NO]/dT */
    J[250] += dqdT;               /* dwdot[N2O]/dT */
    J[251] -= dqdT;               /* dwdot[NH]/dT */

    /*reaction 45: NH + NO <=> N2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = prefactor_units[44] * fwd_A[44]
                * exp(fwd_beta[44] * tc[0] - activation_units[44] * fwd_Ea[44] * invT);
    dlnkfdT = fwd_beta[44] * invT + activation_units[44] * fwd_Ea[44] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[10];
    Kc = exp(-g_RT[0] + g_RT[8] - g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (h_RT[0] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[8] -= q; /* NO */
    wdot[10] += q; /* N2O */
    wdot[11] -= q; /* NH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[10];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[8] -= dqdci;                /* dwdot[NO]/d[H] */
    J[10] += dqdci;               /* dwdot[N2O]/d[H] */
    J[11] -= dqdci;               /* dwdot[NH]/d[H] */
    /* d()/d[NO] */
    dqdci =  + k_f*sc[11];
    J[128] += dqdci;              /* dwdot[H]/d[NO] */
    J[136] -= dqdci;              /* dwdot[NO]/d[NO] */
    J[138] += dqdci;              /* dwdot[N2O]/d[NO] */
    J[139] -= dqdci;              /* dwdot[NH]/d[NO] */
    /* d()/d[N2O] */
    dqdci =  - k_r*sc[0];
    J[160] += dqdci;              /* dwdot[H]/d[N2O] */
    J[168] -= dqdci;              /* dwdot[NO]/d[N2O] */
    J[170] += dqdci;              /* dwdot[N2O]/d[N2O] */
    J[171] -= dqdci;              /* dwdot[NH]/d[N2O] */
    /* d()/d[NH] */
    dqdci =  + k_f*sc[8];
    J[176] += dqdci;              /* dwdot[H]/d[NH] */
    J[184] -= dqdci;              /* dwdot[NO]/d[NH] */
    J[186] += dqdci;              /* dwdot[N2O]/d[NH] */
    J[187] -= dqdci;              /* dwdot[NH]/d[NH] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[248] -= dqdT;               /* dwdot[NO]/dT */
    J[250] += dqdT;               /* dwdot[N2O]/dT */
    J[251] -= dqdT;               /* dwdot[NH]/dT */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = prefactor_units[45] * fwd_A[45]
                * exp(fwd_beta[45] * tc[0] - activation_units[45] * fwd_Ea[45] * invT);
    dlnkfdT = fwd_beta[45] * invT + activation_units[45] * fwd_Ea[45] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[14];
    Kc = exp(-g_RT[2] + g_RT[8] + g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (h_RT[2] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* OH */
    wdot[8] -= q; /* NO */
    wdot[11] -= q; /* NH */
    wdot[14] += q; /* N2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[40] -= dqdci;               /* dwdot[NO]/d[OH] */
    J[43] -= dqdci;               /* dwdot[NH]/d[OH] */
    J[46] += dqdci;               /* dwdot[N2]/d[OH] */
    /* d()/d[NO] */
    dqdci =  + k_f*sc[11];
    J[130] += dqdci;              /* dwdot[OH]/d[NO] */
    J[136] -= dqdci;              /* dwdot[NO]/d[NO] */
    J[139] -= dqdci;              /* dwdot[NH]/d[NO] */
    J[142] += dqdci;              /* dwdot[N2]/d[NO] */
    /* d()/d[NH] */
    dqdci =  + k_f*sc[8];
    J[178] += dqdci;              /* dwdot[OH]/d[NH] */
    J[184] -= dqdci;              /* dwdot[NO]/d[NH] */
    J[187] -= dqdci;              /* dwdot[NH]/d[NH] */
    J[190] += dqdci;              /* dwdot[N2]/d[NH] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[2];
    J[226] += dqdci;              /* dwdot[OH]/d[N2] */
    J[232] -= dqdci;              /* dwdot[NO]/d[N2] */
    J[235] -= dqdci;              /* dwdot[NH]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[248] -= dqdT;               /* dwdot[NO]/dT */
    J[251] -= dqdT;               /* dwdot[NH]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[11];
    k_f = prefactor_units[46] * fwd_A[46]
                * exp(fwd_beta[46] * tc[0] - activation_units[46] * fwd_Ea[46] * invT);
    dlnkfdT = fwd_beta[46] * invT + activation_units[46] * fwd_Ea[46] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[10];
    Kc = exp(-g_RT[2] + g_RT[9] - g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[11]) + (h_RT[2] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* OH */
    wdot[9] -= q; /* NO2 */
    wdot[10] += q; /* N2O */
    wdot[11] -= q; /* NH */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[10];
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[41] -= dqdci;               /* dwdot[NO2]/d[OH] */
    J[42] += dqdci;               /* dwdot[N2O]/d[OH] */
    J[43] -= dqdci;               /* dwdot[NH]/d[OH] */
    /* d()/d[NO2] */
    dqdci =  + k_f*sc[11];
    J[146] += dqdci;              /* dwdot[OH]/d[NO2] */
    J[153] -= dqdci;              /* dwdot[NO2]/d[NO2] */
    J[154] += dqdci;              /* dwdot[N2O]/d[NO2] */
    J[155] -= dqdci;              /* dwdot[NH]/d[NO2] */
    /* d()/d[N2O] */
    dqdci =  - k_r*sc[2];
    J[162] += dqdci;              /* dwdot[OH]/d[N2O] */
    J[169] -= dqdci;              /* dwdot[NO2]/d[N2O] */
    J[170] += dqdci;              /* dwdot[N2O]/d[N2O] */
    J[171] -= dqdci;              /* dwdot[NH]/d[N2O] */
    /* d()/d[NH] */
    dqdci =  + k_f*sc[9];
    J[178] += dqdci;              /* dwdot[OH]/d[NH] */
    J[185] -= dqdci;              /* dwdot[NO2]/d[NH] */
    J[186] += dqdci;              /* dwdot[N2O]/d[NH] */
    J[187] -= dqdci;              /* dwdot[NH]/d[NH] */
    /* d()/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[249] -= dqdT;               /* dwdot[NO2]/dT */
    J[250] += dqdT;               /* dwdot[N2O]/dT */
    J[251] -= dqdT;               /* dwdot[NH]/dT */

    /*reaction 48: N + OH <=> NO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[12];
    k_f = prefactor_units[47] * fwd_A[47]
                * exp(fwd_beta[47] * tc[0] - activation_units[47] * fwd_Ea[47] * invT);
    dlnkfdT = fwd_beta[47] * invT + activation_units[47] * fwd_Ea[47] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[8];
    Kc = exp(-g_RT[0] + g_RT[2] - g_RT[8] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[12]) + (h_RT[0] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[2] -= q; /* OH */
    wdot[8] += q; /* NO */
    wdot[12] -= q; /* N */
    /* d()/d[H] */
    dqdci =  - k_r*sc[8];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[2] -= dqdci;                /* dwdot[OH]/d[H] */
    J[8] += dqdci;                /* dwdot[NO]/d[H] */
    J[12] -= dqdci;               /* dwdot[N]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[32] += dqdci;               /* dwdot[H]/d[OH] */
    J[34] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[40] += dqdci;               /* dwdot[NO]/d[OH] */
    J[44] -= dqdci;               /* dwdot[N]/d[OH] */
    /* d()/d[NO] */
    dqdci =  - k_r*sc[0];
    J[128] += dqdci;              /* dwdot[H]/d[NO] */
    J[130] -= dqdci;              /* dwdot[OH]/d[NO] */
    J[136] += dqdci;              /* dwdot[NO]/d[NO] */
    J[140] -= dqdci;              /* dwdot[N]/d[NO] */
    /* d()/d[N] */
    dqdci =  + k_f*sc[2];
    J[192] += dqdci;              /* dwdot[H]/d[N] */
    J[194] -= dqdci;              /* dwdot[OH]/d[N] */
    J[200] += dqdci;              /* dwdot[NO]/d[N] */
    J[204] -= dqdci;              /* dwdot[N]/d[N] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[242] -= dqdT;               /* dwdot[OH]/dT */
    J[248] += dqdT;               /* dwdot[NO]/dT */
    J[252] -= dqdT;               /* dwdot[N]/dT */

    /*reaction 49: N + O2 <=> NO + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[12];
    k_f = prefactor_units[48] * fwd_A[48]
                * exp(fwd_beta[48] * tc[0] - activation_units[48] * fwd_Ea[48] * invT);
    dlnkfdT = fwd_beta[48] * invT + activation_units[48] * fwd_Ea[48] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[8];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[8] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[12]) + (h_RT[1] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[4] -= q; /* O2 */
    wdot[8] += q; /* NO */
    wdot[12] -= q; /* N */
    /* d()/d[O] */
    dqdci =  - k_r*sc[8];
    J[17] += dqdci;               /* dwdot[O]/d[O] */
    J[20] -= dqdci;               /* dwdot[O2]/d[O] */
    J[24] += dqdci;               /* dwdot[NO]/d[O] */
    J[28] -= dqdci;               /* dwdot[N]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[12];
    J[65] += dqdci;               /* dwdot[O]/d[O2] */
    J[68] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[NO]/d[O2] */
    J[76] -= dqdci;               /* dwdot[N]/d[O2] */
    /* d()/d[NO] */
    dqdci =  - k_r*sc[1];
    J[129] += dqdci;              /* dwdot[O]/d[NO] */
    J[132] -= dqdci;              /* dwdot[O2]/d[NO] */
    J[136] += dqdci;              /* dwdot[NO]/d[NO] */
    J[140] -= dqdci;              /* dwdot[N]/d[NO] */
    /* d()/d[N] */
    dqdci =  + k_f*sc[4];
    J[193] += dqdci;              /* dwdot[O]/d[N] */
    J[196] -= dqdci;              /* dwdot[O2]/d[N] */
    J[200] += dqdci;              /* dwdot[NO]/d[N] */
    J[204] -= dqdci;              /* dwdot[N]/d[N] */
    /* d()/dT */
    J[241] += dqdT;               /* dwdot[O]/dT */
    J[244] -= dqdT;               /* dwdot[O2]/dT */
    J[248] += dqdT;               /* dwdot[NO]/dT */
    J[252] -= dqdT;               /* dwdot[N]/dT */

    /*reaction 50: N + NO <=> N2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = prefactor_units[49] * fwd_A[49]
                * exp(fwd_beta[49] * tc[0] - activation_units[49] * fwd_Ea[49] * invT);
    dlnkfdT = fwd_beta[49] * invT + activation_units[49] * fwd_Ea[49] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[8] + g_RT[12] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[8] -= q; /* NO */
    wdot[12] -= q; /* N */
    wdot[14] += q; /* N2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[14];
    J[17] += dqdci;               /* dwdot[O]/d[O] */
    J[24] -= dqdci;               /* dwdot[NO]/d[O] */
    J[28] -= dqdci;               /* dwdot[N]/d[O] */
    J[30] += dqdci;               /* dwdot[N2]/d[O] */
    /* d()/d[NO] */
    dqdci =  + k_f*sc[12];
    J[129] += dqdci;              /* dwdot[O]/d[NO] */
    J[136] -= dqdci;              /* dwdot[NO]/d[NO] */
    J[140] -= dqdci;              /* dwdot[N]/d[NO] */
    J[142] += dqdci;              /* dwdot[N2]/d[NO] */
    /* d()/d[N] */
    dqdci =  + k_f*sc[8];
    J[193] += dqdci;              /* dwdot[O]/d[N] */
    J[200] -= dqdci;              /* dwdot[NO]/d[N] */
    J[204] -= dqdci;              /* dwdot[N]/d[N] */
    J[206] += dqdci;              /* dwdot[N2]/d[N] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[1];
    J[225] += dqdci;              /* dwdot[O]/d[N2] */
    J[232] -= dqdci;              /* dwdot[NO]/d[N2] */
    J[236] -= dqdci;              /* dwdot[N]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[241] += dqdT;               /* dwdot[O]/dT */
    J[248] -= dqdT;               /* dwdot[NO]/dT */
    J[252] -= dqdT;               /* dwdot[N]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 51: NNH <=> N2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13];
    k_f = prefactor_units[50] * fwd_A[50]
                * exp(fwd_beta[50] * tc[0] - activation_units[50] * fwd_Ea[50] * invT);
    dlnkfdT = fwd_beta[50] * invT + activation_units[50] * fwd_Ea[50] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[14];
    Kc = refC * exp(-g_RT[0] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13]) + (h_RT[0] + h_RT[14]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[13] -= q; /* NNH */
    wdot[14] += q; /* N2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[NNH]/d[H] */
    J[14] += dqdci;               /* dwdot[N2]/d[H] */
    /* d()/d[NNH] */
    dqdci =  + k_f;
    J[208] += dqdci;              /* dwdot[H]/d[NNH] */
    J[221] -= dqdci;              /* dwdot[NNH]/d[NNH] */
    J[222] += dqdci;              /* dwdot[N2]/d[NNH] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[0];
    J[224] += dqdci;              /* dwdot[H]/d[N2] */
    J[237] -= dqdci;              /* dwdot[NNH]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[253] -= dqdT;               /* dwdot[NNH]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[13];
    k_f = prefactor_units[51] * fwd_A[51]
                * exp(fwd_beta[51] * tc[0] - activation_units[51] * fwd_Ea[51] * invT);
    dlnkfdT = fwd_beta[51] * invT + activation_units[51] * fwd_Ea[51] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[14];
    Kc = exp(g_RT[0] - g_RT[3] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[13]) + (h_RT[3] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* H2 */
    wdot[13] -= q; /* NNH */
    wdot[14] += q; /* N2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[13];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[H2]/d[H] */
    J[13] -= dqdci;               /* dwdot[NNH]/d[H] */
    J[14] += dqdci;               /* dwdot[N2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[14];
    J[48] -= dqdci;               /* dwdot[H]/d[H2] */
    J[51] += dqdci;               /* dwdot[H2]/d[H2] */
    J[61] -= dqdci;               /* dwdot[NNH]/d[H2] */
    J[62] += dqdci;               /* dwdot[N2]/d[H2] */
    /* d()/d[NNH] */
    dqdci =  + k_f*sc[0];
    J[208] -= dqdci;              /* dwdot[H]/d[NNH] */
    J[211] += dqdci;              /* dwdot[H2]/d[NNH] */
    J[221] -= dqdci;              /* dwdot[NNH]/d[NNH] */
    J[222] += dqdci;              /* dwdot[N2]/d[NNH] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[3];
    J[224] -= dqdci;              /* dwdot[H]/d[N2] */
    J[227] += dqdci;              /* dwdot[H2]/d[N2] */
    J[237] -= dqdci;              /* dwdot[NNH]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[240] -= dqdT;               /* dwdot[H]/dT */
    J[243] += dqdT;               /* dwdot[H2]/dT */
    J[253] -= dqdT;               /* dwdot[NNH]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 53: NNH + O <=> N2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = prefactor_units[52] * fwd_A[52]
                * exp(fwd_beta[52] * tc[0] - activation_units[52] * fwd_Ea[52] * invT);
    dlnkfdT = fwd_beta[52] * invT + activation_units[52] * fwd_Ea[52] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[10];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[10] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[0] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* O */
    wdot[10] += q; /* N2O */
    wdot[13] -= q; /* NNH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[10];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[O]/d[H] */
    J[10] += dqdci;               /* dwdot[N2O]/d[H] */
    J[13] -= dqdci;               /* dwdot[NNH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[16] += dqdci;               /* dwdot[H]/d[O] */
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[26] += dqdci;               /* dwdot[N2O]/d[O] */
    J[29] -= dqdci;               /* dwdot[NNH]/d[O] */
    /* d()/d[N2O] */
    dqdci =  - k_r*sc[0];
    J[160] += dqdci;              /* dwdot[H]/d[N2O] */
    J[161] -= dqdci;              /* dwdot[O]/d[N2O] */
    J[170] += dqdci;              /* dwdot[N2O]/d[N2O] */
    J[173] -= dqdci;              /* dwdot[NNH]/d[N2O] */
    /* d()/d[NNH] */
    dqdci =  + k_f*sc[1];
    J[208] += dqdci;              /* dwdot[H]/d[NNH] */
    J[209] -= dqdci;              /* dwdot[O]/d[NNH] */
    J[218] += dqdci;              /* dwdot[N2O]/d[NNH] */
    J[221] -= dqdci;              /* dwdot[NNH]/d[NNH] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[250] += dqdT;               /* dwdot[N2O]/dT */
    J[253] -= dqdT;               /* dwdot[NNH]/dT */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = prefactor_units[53] * fwd_A[53]
                * exp(fwd_beta[53] * tc[0] - activation_units[53] * fwd_Ea[53] * invT);
    dlnkfdT = fwd_beta[53] * invT + activation_units[53] * fwd_Ea[53] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[14];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[2] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* OH */
    wdot[13] -= q; /* NNH */
    wdot[14] += q; /* N2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[18] += dqdci;               /* dwdot[OH]/d[O] */
    J[29] -= dqdci;               /* dwdot[NNH]/d[O] */
    J[30] += dqdci;               /* dwdot[N2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[33] -= dqdci;               /* dwdot[O]/d[OH] */
    J[34] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] -= dqdci;               /* dwdot[NNH]/d[OH] */
    J[46] += dqdci;               /* dwdot[N2]/d[OH] */
    /* d()/d[NNH] */
    dqdci =  + k_f*sc[1];
    J[209] -= dqdci;              /* dwdot[O]/d[NNH] */
    J[210] += dqdci;              /* dwdot[OH]/d[NNH] */
    J[221] -= dqdci;              /* dwdot[NNH]/d[NNH] */
    J[222] += dqdci;              /* dwdot[N2]/d[NNH] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[2];
    J[225] -= dqdci;              /* dwdot[O]/d[N2] */
    J[226] += dqdci;              /* dwdot[OH]/d[N2] */
    J[237] -= dqdci;              /* dwdot[NNH]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[242] += dqdT;               /* dwdot[OH]/dT */
    J[253] -= dqdT;               /* dwdot[NNH]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 55: NNH + O <=> NH + NO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = prefactor_units[54] * fwd_A[54]
                * exp(fwd_beta[54] * tc[0] - activation_units[54] * fwd_Ea[54] * invT);
    dlnkfdT = fwd_beta[54] * invT + activation_units[54] * fwd_Ea[54] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[11];
    Kc = exp(g_RT[1] - g_RT[8] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[8] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[8] += q; /* NO */
    wdot[11] += q; /* NH */
    wdot[13] -= q; /* NNH */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[17] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += dqdci;               /* dwdot[NO]/d[O] */
    J[27] += dqdci;               /* dwdot[NH]/d[O] */
    J[29] -= dqdci;               /* dwdot[NNH]/d[O] */
    /* d()/d[NO] */
    dqdci =  - k_r*sc[11];
    J[129] -= dqdci;              /* dwdot[O]/d[NO] */
    J[136] += dqdci;              /* dwdot[NO]/d[NO] */
    J[139] += dqdci;              /* dwdot[NH]/d[NO] */
    J[141] -= dqdci;              /* dwdot[NNH]/d[NO] */
    /* d()/d[NH] */
    dqdci =  - k_r*sc[8];
    J[177] -= dqdci;              /* dwdot[O]/d[NH] */
    J[184] += dqdci;              /* dwdot[NO]/d[NH] */
    J[187] += dqdci;              /* dwdot[NH]/d[NH] */
    J[189] -= dqdci;              /* dwdot[NNH]/d[NH] */
    /* d()/d[NNH] */
    dqdci =  + k_f*sc[1];
    J[209] -= dqdci;              /* dwdot[O]/d[NNH] */
    J[216] += dqdci;              /* dwdot[NO]/d[NNH] */
    J[219] += dqdci;              /* dwdot[NH]/d[NNH] */
    J[221] -= dqdci;              /* dwdot[NNH]/d[NNH] */
    /* d()/dT */
    J[241] -= dqdT;               /* dwdot[O]/dT */
    J[248] += dqdT;               /* dwdot[NO]/dT */
    J[251] += dqdT;               /* dwdot[NH]/dT */
    J[253] -= dqdT;               /* dwdot[NNH]/dT */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = prefactor_units[55] * fwd_A[55]
                * exp(fwd_beta[55] * tc[0] - activation_units[55] * fwd_Ea[55] * invT);
    dlnkfdT = fwd_beta[55] * invT + activation_units[55] * fwd_Ea[55] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[14];
    Kc = exp(g_RT[2] - g_RT[6] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[6] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[13] -= q; /* NNH */
    wdot[14] += q; /* N2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[34] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[38] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[45] -= dqdci;               /* dwdot[NNH]/d[OH] */
    J[46] += dqdci;               /* dwdot[N2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[14];
    J[98] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[109] -= dqdci;              /* dwdot[NNH]/d[H2O] */
    J[110] += dqdci;              /* dwdot[N2]/d[H2O] */
    /* d()/d[NNH] */
    dqdci =  + k_f*sc[2];
    J[210] -= dqdci;              /* dwdot[OH]/d[NNH] */
    J[214] += dqdci;              /* dwdot[H2O]/d[NNH] */
    J[221] -= dqdci;              /* dwdot[NNH]/d[NNH] */
    J[222] += dqdci;              /* dwdot[N2]/d[NNH] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[6];
    J[226] -= dqdci;              /* dwdot[OH]/d[N2] */
    J[230] += dqdci;              /* dwdot[H2O]/d[N2] */
    J[237] -= dqdci;              /* dwdot[NNH]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[242] -= dqdT;               /* dwdot[OH]/dT */
    J[246] += dqdT;               /* dwdot[H2O]/dT */
    J[253] -= dqdT;               /* dwdot[NNH]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = prefactor_units[56] * fwd_A[56]
                * exp(fwd_beta[56] * tc[0] - activation_units[56] * fwd_Ea[56] * invT);
    dlnkfdT = fwd_beta[56] * invT + activation_units[56] * fwd_Ea[56] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[14];
    Kc = exp(g_RT[4] - g_RT[5] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[5] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O2 */
    wdot[5] += q; /* HO2 */
    wdot[13] -= q; /* NNH */
    wdot[14] += q; /* N2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[68] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[69] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[77] -= dqdci;               /* dwdot[NNH]/d[O2] */
    J[78] += dqdci;               /* dwdot[N2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[14];
    J[84] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[85] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[93] -= dqdci;               /* dwdot[NNH]/d[HO2] */
    J[94] += dqdci;               /* dwdot[N2]/d[HO2] */
    /* d()/d[NNH] */
    dqdci =  + k_f*sc[4];
    J[212] -= dqdci;              /* dwdot[O2]/d[NNH] */
    J[213] += dqdci;              /* dwdot[HO2]/d[NNH] */
    J[221] -= dqdci;              /* dwdot[NNH]/d[NNH] */
    J[222] += dqdci;              /* dwdot[N2]/d[NNH] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[5];
    J[228] -= dqdci;              /* dwdot[O2]/d[N2] */
    J[229] += dqdci;              /* dwdot[HO2]/d[N2] */
    J[237] -= dqdci;              /* dwdot[NNH]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[244] -= dqdT;               /* dwdot[O2]/dT */
    J[245] += dqdT;               /* dwdot[HO2]/dT */
    J[253] -= dqdT;               /* dwdot[NNH]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = prefactor_units[57] * fwd_A[57]
                * exp(fwd_beta[57] * tc[0] - activation_units[57] * fwd_Ea[57] * invT);
    dlnkfdT = fwd_beta[57] * invT + activation_units[57] * fwd_Ea[57] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[4]*sc[14];
    Kc = refC * exp(-g_RT[0] + g_RT[4] - g_RT[4] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[0] + h_RT[4] + h_RT[14]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[13] -= q; /* NNH */
    wdot[14] += q; /* N2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4]*sc[14];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[NNH]/d[H] */
    J[14] += dqdci;               /* dwdot[N2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13] - k_r*sc[0]*sc[14];
    J[64] += dqdci;               /* dwdot[H]/d[O2] */
    J[77] -= dqdci;               /* dwdot[NNH]/d[O2] */
    J[78] += dqdci;               /* dwdot[N2]/d[O2] */
    /* d()/d[NNH] */
    dqdci =  + k_f*sc[4];
    J[208] += dqdci;              /* dwdot[H]/d[NNH] */
    J[221] -= dqdci;              /* dwdot[NNH]/d[NNH] */
    J[222] += dqdci;              /* dwdot[N2]/d[NNH] */
    /* d()/d[N2] */
    dqdci =  - k_r*sc[0]*sc[4];
    J[224] += dqdci;              /* dwdot[H]/d[N2] */
    J[237] -= dqdci;              /* dwdot[NNH]/d[N2] */
    J[238] += dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[240] += dqdT;               /* dwdot[H]/dT */
    J[253] -= dqdT;               /* dwdot[NNH]/dT */
    J[254] += dqdT;               /* dwdot[N2]/dT */

    double c_R[15], dcRdT[15], e_RT[15];
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
    for (int k = 0; k < 15; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[240+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 15; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 15; ++m) {
            dehmixdc += eh_RT[m]*J[k*16+m];
        }
        J[k*16+15] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[255] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void dcvpRdT(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: O */
        species[1] =
            -3.27931884e-03
            +1.32861279e-05 * tc[1]
            -1.83841987e-08 * tc[2]
            +8.45063884e-12 * tc[3];
        /*species 2: OH */
        species[2] =
            -2.40106655e-03
            +9.23328066e-06 * tc[1]
            -1.16374892e-08 * tc[2]
            +5.45278008e-12 * tc[3];
        /*species 3: H2 */
        species[3] =
            +7.98042480e-03
            -3.89558340e-05 * tc[1]
            +6.04709010e-08 * tc[2]
            -2.95041156e-11 * tc[3];
        /*species 4: O2 */
        species[4] =
            -2.99673415e-03
            +1.96946040e-05 * tc[1]
            -2.90438852e-08 * tc[2]
            +1.29749134e-11 * tc[3];
        /*species 5: HO2 */
        species[5] =
            -4.74902010e-03
            +4.23159060e-05 * tc[1]
            -7.28278830e-08 * tc[2]
            +3.71682680e-11 * tc[3];
        /*species 6: H2O */
        species[6] =
            -2.03640170e-03
            +1.30406832e-05 * tc[1]
            -1.64637807e-08 * tc[2]
            +7.08787200e-12 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            -8.47390622e-04
            +3.52808646e-05 * tc[1]
            -6.80288832e-08 * tc[2]
            +3.63580063e-11 * tc[3];
        /*species 8: NO */
        species[8] =
            -4.63988124e-03
            +2.20886098e-05 * tc[1]
            -2.80216652e-08 * tc[2]
            +1.12221950e-11 * tc[3];
        /*species 9: NO2 */
        species[9] =
            -1.58542900e-03
            +3.33156240e-05 * tc[1]
            -6.14262780e-08 * tc[2]
            +3.13402256e-11 * tc[3];
        /*species 10: N2O */
        species[10] =
            +9.49219300e-03
            -1.95855500e-05 * tc[1]
            +1.87915350e-08 * tc[2]
            -7.60730400e-12 * tc[3];
        /*species 11: NH */
        species[11] =
            +3.11791970e-04
            -2.97809680e-06 * tc[1]
            +7.44493260e-09 * tc[2]
            -4.14278680e-12 * tc[3];
        /*species 12: N */
        species[12] =
            -2.18001800e-05
            +1.08410580e-07 * tc[1]
            -1.69426800e-10 * tc[2]
            +8.39961600e-14 * tc[3];
        /*species 13: NNH */
        species[13] =
            +2.05358700e-03
            +1.43408200e-06 * tc[1]
            +1.47640440e-09 * tc[2]
            -3.86846800e-12 * tc[3];
        /*species 14: N2 */
        species[14] =
            -1.23660988e-04
            -1.00599887e-06 * tc[1]
            +7.30591836e-09 * tc[2]
            -5.63524940e-12 * tc[3];
    } else {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: O */
        species[1] =
            -2.73162486e-05
            -8.38059040e-09 * tc[1]
            +1.48644553e-11 * tc[2]
            -1.91821478e-15 * tc[3];
        /*species 2: OH */
        species[2] =
            +1.10741289e-03
            -5.88000418e-07 * tc[1]
            +1.26209619e-10 * tc[2]
            -9.69159560e-15 * tc[3];
        /*species 3: H2 */
        species[3] =
            +8.26598020e-04
            -2.92801140e-07 * tc[1]
            +4.62295530e-11 * tc[2]
            -2.75518460e-15 * tc[3];
        /*species 4: O2 */
        species[4] =
            +6.56365523e-04
            -2.82298970e-07 * tc[1]
            +6.17392974e-11 * tc[2]
            -5.19652992e-15 * tc[3];
        /*species 5: HO2 */
        species[5] =
            +1.88120980e-03
            -6.92585940e-07 * tc[1]
            +5.84055480e-11 * tc[2]
            +7.04366120e-16 * tc[3];
        /*species 6: H2O */
        species[6] =
            +2.97318160e-03
            -1.54753778e-06 * tc[1]
            +2.83300542e-10 * tc[2]
            -1.70759964e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.05326003e-03
            -2.59689460e-06 * tc[1]
            +5.94634200e-10 * tc[2]
            -4.55875168e-14 * tc[3];
        /*species 8: NO */
        species[8] =
            +1.19101135e-03
            -8.58245292e-07 * tc[1]
            +2.08344439e-10 * tc[2]
            -1.61318272e-14 * tc[3];
        /*species 9: NO2 */
        species[9] =
            +2.17239550e-03
            -1.65613818e-06 * tc[1]
            +4.72425300e-10 * tc[2]
            -4.20435800e-14 * tc[3];
        /*species 10: N2O */
        species[10] =
            +2.87371400e-03
            -2.39499200e-06 * tc[1]
            +6.75165600e-10 * tc[2]
            -6.30134800e-14 * tc[3];
        /*species 11: NH */
        species[11] =
            +1.32984290e-03
            -8.49560940e-07 * tc[1]
            +2.35045512e-10 * tc[2]
            -2.20177880e-14 * tc[3];
        /*species 12: N */
        species[12] =
            +1.06614600e-04
            -1.49306740e-07 * tc[1]
            +5.63895600e-11 * tc[2]
            -4.10393600e-15 * tc[3];
        /*species 13: NNH */
        species[13] =
            +1.61438800e-03
            -3.26578800e-07 * tc[1]
            -2.56795380e-10 * tc[2]
            +6.45916400e-14 * tc[3];
        /*species 14: N2 */
        species[14] =
            +1.39690040e-03
            -9.85263206e-07 * tc[1]
            +2.35803059e-10 * tc[2]
            -1.84302082e-14 * tc[3];
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

    double q_f[58], q_r[58];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 58; ++i) {
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
    kc[0] = 1.0 / (refC) * exp((g_RT[0] + g_RT[4]) - (g_RT[5]));

    /*reaction 2: H + O2 (+N2) <=> HO2 (+N2) */
    kc[1] = 1.0 / (refC) * exp((g_RT[0] + g_RT[4]) - (g_RT[5]));

    /*reaction 3: H2O2 (+M) <=> OH + OH (+M) */
    kc[2] = refC * exp((g_RT[7]) - (g_RT[2] + g_RT[2]));

    /*reaction 4: NO + O (+M) <=> NO2 (+M) */
    kc[3] = 1.0 / (refC) * exp((g_RT[8] + g_RT[1]) - (g_RT[9]));

    /*reaction 5: N2O (+M) <=> N2 + O (+M) */
    kc[4] = refC * exp((g_RT[10]) - (g_RT[14] + g_RT[1]));

    /*reaction 6: H + H + M <=> H2 + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0]) - (g_RT[3]));

    /*reaction 7: H + O + M <=> OH + M */
    kc[6] = 1.0 / (refC) * exp((g_RT[0] + g_RT[1]) - (g_RT[2]));

    /*reaction 8: O + O + M <=> O2 + M */
    kc[7] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[4]));

    /*reaction 9: OH + H + M <=> H2O + M */
    kc[8] = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[6]));

    /*reaction 10: H + O2 <=> O + OH */
    kc[9] = exp((g_RT[0] + g_RT[4]) - (g_RT[1] + g_RT[2]));

    /*reaction 11: H + H + N2 <=> H2 + N2 */
    kc[10] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[14]) - (g_RT[3] + g_RT[14]));

    /*reaction 12: H + H + H2 <=> H2 + H2 */
    kc[11] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[3]) - (g_RT[3] + g_RT[3]));

    /*reaction 13: H + H + H2O <=> H2 + H2O */
    kc[12] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[6]) - (g_RT[3] + g_RT[6]));

    /*reaction 14: O + H2 <=> OH + H */
    kc[13] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[0]));

    /*reaction 15: O + H2 <=> OH + H */
    kc[14] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[0]));

    /*reaction 16: OH + OH <=> O + H2O */
    kc[15] = exp((g_RT[2] + g_RT[2]) - (g_RT[1] + g_RT[6]));

    /*reaction 17: OH + H2 <=> H + H2O */
    kc[16] = exp((g_RT[2] + g_RT[3]) - (g_RT[0] + g_RT[6]));

    /*reaction 18: H2 + O2 <=> HO2 + H */
    kc[17] = exp((g_RT[3] + g_RT[4]) - (g_RT[5] + g_RT[0]));

    /*reaction 19: HO2 + H <=> OH + OH */
    kc[18] = exp((g_RT[5] + g_RT[0]) - (g_RT[2] + g_RT[2]));

    /*reaction 20: HO2 + H <=> H2O + O */
    kc[19] = exp((g_RT[5] + g_RT[0]) - (g_RT[6] + g_RT[1]));

    /*reaction 21: HO2 + O <=> OH + O2 */
    kc[20] = exp((g_RT[5] + g_RT[1]) - (g_RT[2] + g_RT[4]));

    /*reaction 22: HO2 + OH <=> H2O + O2 */
    kc[21] = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    kc[22] = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    kc[23] = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    kc[24] = exp((g_RT[5] + g_RT[5]) - (g_RT[7] + g_RT[4]));

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    kc[25] = exp((g_RT[5] + g_RT[5]) - (g_RT[7] + g_RT[4]));

    /*reaction 27: H2O2 + H <=> H2O + OH */
    kc[26] = exp((g_RT[7] + g_RT[0]) - (g_RT[6] + g_RT[2]));

    /*reaction 28: H2O2 + H <=> HO2 + H2 */
    kc[27] = exp((g_RT[7] + g_RT[0]) - (g_RT[5] + g_RT[3]));

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    kc[28] = exp((g_RT[7] + g_RT[1]) - (g_RT[5] + g_RT[2]));

    /*reaction 30: H2O2 + OH <=> H2O + HO2 */
    kc[29] = exp((g_RT[7] + g_RT[2]) - (g_RT[6] + g_RT[5]));

    /*reaction 31: H2O2 + OH <=> H2O + HO2 */
    kc[30] = exp((g_RT[7] + g_RT[2]) - (g_RT[6] + g_RT[5]));

    /*reaction 32: NO + HO2 <=> NO2 + OH */
    kc[31] = exp((g_RT[8] + g_RT[5]) - (g_RT[9] + g_RT[2]));

    /*reaction 33: NO2 + H <=> NO + OH */
    kc[32] = exp((g_RT[9] + g_RT[0]) - (g_RT[8] + g_RT[2]));

    /*reaction 34: NO2 + O <=> NO + O2 */
    kc[33] = exp((g_RT[9] + g_RT[1]) - (g_RT[8] + g_RT[4]));

    /*reaction 35: NO2 + NO2 <=> NO + NO + O2 */
    kc[34] = refC * exp((g_RT[9] + g_RT[9]) - (g_RT[8] + g_RT[8] + g_RT[4]));

    /*reaction 36: N2O + H <=> N2 + OH */
    kc[35] = exp((g_RT[10] + g_RT[0]) - (g_RT[14] + g_RT[2]));

    /*reaction 37: N2O + H <=> N2 + OH */
    kc[36] = exp((g_RT[10] + g_RT[0]) - (g_RT[14] + g_RT[2]));

    /*reaction 38: N2O + O <=> NO + NO */
    kc[37] = exp((g_RT[10] + g_RT[1]) - (g_RT[8] + g_RT[8]));

    /*reaction 39: N2O + O <=> N2 + O2 */
    kc[38] = exp((g_RT[10] + g_RT[1]) - (g_RT[14] + g_RT[4]));

    /*reaction 40: NH + H <=> N + H2 */
    kc[39] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[3]));

    /*reaction 41: NH + O <=> NO + H */
    kc[40] = exp((g_RT[11] + g_RT[1]) - (g_RT[8] + g_RT[0]));

    /*reaction 42: NH + OH <=> N + H2O */
    kc[41] = exp((g_RT[11] + g_RT[2]) - (g_RT[12] + g_RT[6]));

    /*reaction 43: NH + O2 <=> NO + OH */
    kc[42] = exp((g_RT[11] + g_RT[4]) - (g_RT[8] + g_RT[2]));

    /*reaction 44: NH + NO <=> N2O + H */
    kc[43] = exp((g_RT[11] + g_RT[8]) - (g_RT[10] + g_RT[0]));

    /*reaction 45: NH + NO <=> N2O + H */
    kc[44] = exp((g_RT[11] + g_RT[8]) - (g_RT[10] + g_RT[0]));

    /*reaction 46: NH + NO <=> N2 + OH */
    kc[45] = exp((g_RT[11] + g_RT[8]) - (g_RT[14] + g_RT[2]));

    /*reaction 47: NH + NO2 <=> N2O + OH */
    kc[46] = exp((g_RT[11] + g_RT[9]) - (g_RT[10] + g_RT[2]));

    /*reaction 48: N + OH <=> NO + H */
    kc[47] = exp((g_RT[12] + g_RT[2]) - (g_RT[8] + g_RT[0]));

    /*reaction 49: N + O2 <=> NO + O */
    kc[48] = exp((g_RT[12] + g_RT[4]) - (g_RT[8] + g_RT[1]));

    /*reaction 50: N + NO <=> N2 + O */
    kc[49] = exp((g_RT[12] + g_RT[8]) - (g_RT[14] + g_RT[1]));

    /*reaction 51: NNH <=> N2 + H */
    kc[50] = refC * exp((g_RT[13]) - (g_RT[14] + g_RT[0]));

    /*reaction 52: NNH + H <=> N2 + H2 */
    kc[51] = exp((g_RT[13] + g_RT[0]) - (g_RT[14] + g_RT[3]));

    /*reaction 53: NNH + O <=> N2O + H */
    kc[52] = exp((g_RT[13] + g_RT[1]) - (g_RT[10] + g_RT[0]));

    /*reaction 54: NNH + O <=> N2 + OH */
    kc[53] = exp((g_RT[13] + g_RT[1]) - (g_RT[14] + g_RT[2]));

    /*reaction 55: NNH + O <=> NH + NO */
    kc[54] = exp((g_RT[13] + g_RT[1]) - (g_RT[11] + g_RT[8]));

    /*reaction 56: NNH + OH <=> N2 + H2O */
    kc[55] = exp((g_RT[13] + g_RT[2]) - (g_RT[14] + g_RT[6]));

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    kc[56] = exp((g_RT[13] + g_RT[4]) - (g_RT[14] + g_RT[5]));

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    kc[57] = refC * exp((g_RT[13] + g_RT[4]) - (g_RT[14] + g_RT[0] + g_RT[4]));

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
        /*species 0: H */
        species[0] =
            +2.547366000000000e+04 * invT
            +2.946682850000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.371652480000000e+03 * invT
            +4.095798299000000e+00
            -3.991984240000000e+00 * tc[0]
            +1.200533275000000e-03 * tc[1]
            -7.694400550000000e-07 * tc[2]
            +3.232635883333333e-10 * tc[3]
            -6.815975100000000e-14 * tc[4];
        /*species 3: H2 */
        species[3] =
            -9.179241300000000e+02 * invT
            +1.661300720000000e+00
            -2.344302900000000e+00 * tc[0]
            -3.990212400000000e-03 * tc[1]
            +3.246319500000000e-06 * tc[2]
            -1.679747250000000e-09 * tc[3]
            +3.688014450000000e-13 * tc[4];
        /*species 4: O2 */
        species[4] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367075000000e-03 * tc[1]
            -1.641217000000000e-06 * tc[2]
            +8.067745900000000e-10 * tc[3]
            -1.621864180000000e-13 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +2.631909830000000e+02 * invT
            +5.859106000000001e-01
            -4.301788000000000e+00 * tc[0]
            +2.374510050000000e-03 * tc[1]
            -3.526325500000000e-06 * tc[2]
            +2.022996750000000e-09 * tc[3]
            -4.646033500000000e-13 * tc[4];
        /*species 6: H2O */
        species[6] =
            -3.029372600000000e+04 * invT
            +5.047644210000000e+00
            -4.198635200000000e+00 * tc[0]
            +1.018200850000000e-03 * tc[1]
            -1.086723600000000e-06 * tc[2]
            +4.573272416666666e-10 * tc[3]
            -8.859840000000001e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.768436010000000e+04 * invT
            +1.041419330000000e+00
            -4.315151490000000e+00 * tc[0]
            +4.236953110000000e-04 * tc[1]
            -2.940072050000000e-06 * tc[2]
            +1.889691200000000e-09 * tc[3]
            -4.544750790000000e-13 * tc[4];
        /*species 8: NO */
        species[8] =
            +9.818237859999999e+03 * invT
            +1.937989440000000e+00
            -4.218598960000000e+00 * tc[0]
            +2.319940620000000e-03 * tc[1]
            -1.840717483333334e-06 * tc[2]
            +7.783795891666667e-10 * tc[3]
            -1.402774370000000e-13 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +2.874097570000000e+03 * invT
            -2.367960700000000e+00
            -3.944031200000000e+00 * tc[0]
            +7.927145000000000e-04 * tc[1]
            -2.776302000000000e-06 * tc[2]
            +1.706285500000000e-09 * tc[3]
            -3.917528200000000e-13 * tc[4];
        /*species 10: N2O */
        species[10] =
            +8.765100000000000e+03 * invT
            -6.968164000000000e+00
            -2.543058000000000e+00 * tc[0]
            -4.746096500000000e-03 * tc[1]
            +1.632129166666667e-06 * tc[2]
            -5.219870833333333e-10 * tc[3]
            +9.509129999999999e-14 * tc[4];
        /*species 11: NH */
        species[11] =
            +4.189429400000000e+04 * invT
            +1.644580700000000e+00
            -3.492908400000000e+00 * tc[0]
            -1.558959850000000e-04 * tc[1]
            +2.481747333333333e-07 * tc[2]
            -2.068036833333333e-10 * tc[3]
            +5.178483500000000e-14 * tc[4];
        /*species 12: N */
        species[12] =
            +5.609890000000000e+04 * invT
            -1.664495000000000e+00
            -2.503071000000000e+00 * tc[0]
            +1.090009000000000e-05 * tc[1]
            -9.034215000000001e-09 * tc[2]
            +4.706300000000000e-12 * tc[3]
            -1.049952000000000e-15 * tc[4];
        /*species 13: NNH */
        species[13] =
            +2.833347000000000e+04 * invT
            -2.890493000000000e+00
            -3.501344000000000e+00 * tc[0]
            -1.026793500000000e-03 * tc[1]
            -1.195068333333333e-07 * tc[2]
            -4.101123333333334e-11 * tc[3]
            +4.835585000000000e-14 * tc[4];
        /*species 14: N2 */
        species[14] =
            -1.046976280000000e+03 * invT
            +5.635349000000001e-01
            -3.531005280000000e+00 * tc[0]
            +6.183049400000000e-05 * tc[1]
            +8.383323883333334e-08 * tc[2]
            -2.029421766666667e-10 * tc[3]
            +7.044061750000001e-14 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.547366000000000e+04 * invT
            +2.946682850000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.922601200000000e+04 * invT
            -2.378657600000000e+00
            -2.543636970000000e+00 * tc[0]
            +1.365812430000000e-05 * tc[1]
            +6.983825333333333e-10 * tc[2]
            -4.129015375000000e-13 * tc[3]
            +2.397768470000000e-17 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.700562200000000e+03 * invT
            -3.006600610000000e+00
            -2.838530330000000e+00 * tc[0]
            -5.537064450000000e-04 * tc[1]
            +4.900003483333333e-08 * tc[2]
            -3.505822741666666e-12 * tc[3]
            +1.211449450000000e-16 * tc[4];
        /*species 3: H2 */
        species[3] =
            -8.130558200000000e+02 * invT
            +3.957146900000000e+00
            -2.932830500000000e+00 * tc[0]
            -4.132990100000000e-04 * tc[1]
            +2.440009500000000e-08 * tc[2]
            -1.284154250000000e-12 * tc[3]
            +3.443980750000000e-17 * tc[4];
        /*species 4: O2 */
        species[4] =
            -1.215977250000000e+03 * invT
            +2.455989900000000e-01
            -3.660960830000000e+00 * tc[0]
            -3.281827615000000e-04 * tc[1]
            +2.352491416666667e-08 * tc[2]
            -1.714980483333333e-12 * tc[3]
            +6.495662400000000e-17 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +3.020107360000000e+01 * invT
            +1.215292100000000e+00
            -4.172265900000000e+00 * tc[0]
            -9.406049000000000e-04 * tc[1]
            +5.771549500000000e-08 * tc[2]
            -1.622376333333333e-12 * tc[3]
            -8.804576500000000e-18 * tc[4];
        /*species 6: H2O */
        species[6] =
            -2.988589400000000e+04 * invT
            -4.205511100000001e+00
            -2.677038900000000e+00 * tc[0]
            -1.486590800000000e-03 * tc[1]
            +1.289614816666667e-07 * tc[2]
            -7.869459500000001e-12 * tc[3]
            +2.134499550000000e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.798479390000000e+04 * invT
            +3.914803390000000e+00
            -4.579773050000000e+00 * tc[0]
            -2.026630015000000e-03 * tc[1]
            +2.164078833333333e-07 * tc[2]
            -1.651761666666667e-11 * tc[3]
            +5.698439600000000e-16 * tc[4];
        /*species 8: NO */
        species[8] =
            +9.894569540000000e+03 * invT
            -3.108292350000000e+00
            -3.260712340000000e+00 * tc[0]
            -5.955056750000000e-04 * tc[1]
            +7.152044099999999e-08 * tc[2]
            -5.787345524999999e-12 * tc[3]
            +2.016478405000000e-16 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +2.293977770000000e+03 * invT
            +5.002170951000000e+00
            -4.884754000000000e+00 * tc[0]
            -1.086197750000000e-03 * tc[1]
            +1.380115150000000e-07 * tc[2]
            -1.312292500000000e-11 * tc[3]
            +5.255447500000000e-16 * tc[4];
        /*species 10: N2O */
        species[10] =
            +8.165811000000000e+03 * invT
            +6.376227000000000e+00
            -4.718977000000000e+00 * tc[0]
            -1.436857000000000e-03 * tc[1]
            +1.995826666666667e-07 * tc[2]
            -1.875460000000000e-11 * tc[3]
            +7.876685000000001e-16 * tc[4];
        /*species 11: NH */
        species[11] =
            +4.213451400000000e+04 * invT
            -2.957086900000000e+00
            -2.783692900000000e+00 * tc[0]
            -6.649214500000000e-04 * tc[1]
            +7.079674500000000e-08 * tc[2]
            -6.529041999999999e-12 * tc[3]
            +2.752223500000000e-16 * tc[4];
        /*species 12: N */
        species[12] =
            +5.611604000000000e+04 * invT
            -1.998490000000000e+00
            -2.450268000000000e+00 * tc[0]
            -5.330730000000000e-05 * tc[1]
            +1.244222833333333e-08 * tc[2]
            -1.566376666666667e-12 * tc[3]
            +5.129920000000000e-17 * tc[4];
        /*species 13: NNH */
        species[13] =
            +2.788029000000000e+04 * invT
            +3.511053200000000e+00
            -4.415342000000000e+00 * tc[0]
            -8.071940000000000e-04 * tc[1]
            +2.721490000000000e-08 * tc[2]
            +7.133205000000000e-12 * tc[3]
            -8.073955000000000e-16 * tc[4];
        /*species 14: N2 */
        species[14] =
            -9.239486879999999e+02 * invT
            -2.919311250000000e+00
            -2.952576370000000e+00 * tc[0]
            -6.984502000000001e-04 * tc[1]
            +8.210526716666667e-08 * tc[2]
            -6.550084958333333e-12 * tc[3]
            +2.303776020000000e-16 * tc[4];
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
        /*species 0: H */
        species[0] =
            +2.54736600e+04 * invT
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.91222592e+04 * invT
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.37165248e+03 * invT
            +3.09579830e+00
            -3.99198424e+00 * tc[0]
            +1.20053327e-03 * tc[1]
            -7.69440055e-07 * tc[2]
            +3.23263588e-10 * tc[3]
            -6.81597510e-14 * tc[4];
        /*species 3: H2 */
        species[3] =
            -9.17924130e+02 * invT
            +6.61300720e-01
            -2.34430290e+00 * tc[0]
            -3.99021240e-03 * tc[1]
            +3.24631950e-06 * tc[2]
            -1.67974725e-09 * tc[3]
            +3.68801445e-13 * tc[4];
        /*species 4: O2 */
        species[4] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836707e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774590e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +2.63190983e+02 * invT
            -4.14089400e-01
            -4.30178800e+00 * tc[0]
            +2.37451005e-03 * tc[1]
            -3.52632550e-06 * tc[2]
            +2.02299675e-09 * tc[3]
            -4.64603350e-13 * tc[4];
        /*species 6: H2O */
        species[6] =
            -3.02937260e+04 * invT
            +4.04764421e+00
            -4.19863520e+00 * tc[0]
            +1.01820085e-03 * tc[1]
            -1.08672360e-06 * tc[2]
            +4.57327242e-10 * tc[3]
            -8.85984000e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76843601e+04 * invT
            +4.14193300e-02
            -4.31515149e+00 * tc[0]
            +4.23695311e-04 * tc[1]
            -2.94007205e-06 * tc[2]
            +1.88969120e-09 * tc[3]
            -4.54475079e-13 * tc[4];
        /*species 8: NO */
        species[8] =
            +9.81823786e+03 * invT
            +9.37989440e-01
            -4.21859896e+00 * tc[0]
            +2.31994062e-03 * tc[1]
            -1.84071748e-06 * tc[2]
            +7.78379589e-10 * tc[3]
            -1.40277437e-13 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +2.87409757e+03 * invT
            -3.36796070e+00
            -3.94403120e+00 * tc[0]
            +7.92714500e-04 * tc[1]
            -2.77630200e-06 * tc[2]
            +1.70628550e-09 * tc[3]
            -3.91752820e-13 * tc[4];
        /*species 10: N2O */
        species[10] =
            +8.76510000e+03 * invT
            -7.96816400e+00
            -2.54305800e+00 * tc[0]
            -4.74609650e-03 * tc[1]
            +1.63212917e-06 * tc[2]
            -5.21987083e-10 * tc[3]
            +9.50913000e-14 * tc[4];
        /*species 11: NH */
        species[11] =
            +4.18942940e+04 * invT
            +6.44580700e-01
            -3.49290840e+00 * tc[0]
            -1.55895985e-04 * tc[1]
            +2.48174733e-07 * tc[2]
            -2.06803683e-10 * tc[3]
            +5.17848350e-14 * tc[4];
        /*species 12: N */
        species[12] =
            +5.60989000e+04 * invT
            -2.66449500e+00
            -2.50307100e+00 * tc[0]
            +1.09000900e-05 * tc[1]
            -9.03421500e-09 * tc[2]
            +4.70630000e-12 * tc[3]
            -1.04995200e-15 * tc[4];
        /*species 13: NNH */
        species[13] =
            +2.83334700e+04 * invT
            -3.89049300e+00
            -3.50134400e+00 * tc[0]
            -1.02679350e-03 * tc[1]
            -1.19506833e-07 * tc[2]
            -4.10112333e-11 * tc[3]
            +4.83558500e-14 * tc[4];
        /*species 14: N2 */
        species[14] =
            -1.04697628e+03 * invT
            -4.36465100e-01
            -3.53100528e+00 * tc[0]
            +6.18304940e-05 * tc[1]
            +8.38332388e-08 * tc[2]
            -2.02942177e-10 * tc[3]
            +7.04406175e-14 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.54736600e+04 * invT
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.92260120e+04 * invT
            -3.37865760e+00
            -2.54363697e+00 * tc[0]
            +1.36581243e-05 * tc[1]
            +6.98382533e-10 * tc[2]
            -4.12901538e-13 * tc[3]
            +2.39776847e-17 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.70056220e+03 * invT
            -4.00660061e+00
            -2.83853033e+00 * tc[0]
            -5.53706445e-04 * tc[1]
            +4.90000348e-08 * tc[2]
            -3.50582274e-12 * tc[3]
            +1.21144945e-16 * tc[4];
        /*species 3: H2 */
        species[3] =
            -8.13055820e+02 * invT
            +2.95714690e+00
            -2.93283050e+00 * tc[0]
            -4.13299010e-04 * tc[1]
            +2.44000950e-08 * tc[2]
            -1.28415425e-12 * tc[3]
            +3.44398075e-17 * tc[4];
        /*species 4: O2 */
        species[4] =
            -1.21597725e+03 * invT
            -7.54401010e-01
            -3.66096083e+00 * tc[0]
            -3.28182761e-04 * tc[1]
            +2.35249142e-08 * tc[2]
            -1.71498048e-12 * tc[3]
            +6.49566240e-17 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +3.02010736e+01 * invT
            +2.15292100e-01
            -4.17226590e+00 * tc[0]
            -9.40604900e-04 * tc[1]
            +5.77154950e-08 * tc[2]
            -1.62237633e-12 * tc[3]
            -8.80457650e-18 * tc[4];
        /*species 6: H2O */
        species[6] =
            -2.98858940e+04 * invT
            -5.20551110e+00
            -2.67703890e+00 * tc[0]
            -1.48659080e-03 * tc[1]
            +1.28961482e-07 * tc[2]
            -7.86945950e-12 * tc[3]
            +2.13449955e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.79847939e+04 * invT
            +2.91480339e+00
            -4.57977305e+00 * tc[0]
            -2.02663002e-03 * tc[1]
            +2.16407883e-07 * tc[2]
            -1.65176167e-11 * tc[3]
            +5.69843960e-16 * tc[4];
        /*species 8: NO */
        species[8] =
            +9.89456954e+03 * invT
            -4.10829235e+00
            -3.26071234e+00 * tc[0]
            -5.95505675e-04 * tc[1]
            +7.15204410e-08 * tc[2]
            -5.78734552e-12 * tc[3]
            +2.01647840e-16 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +2.29397777e+03 * invT
            +4.00217095e+00
            -4.88475400e+00 * tc[0]
            -1.08619775e-03 * tc[1]
            +1.38011515e-07 * tc[2]
            -1.31229250e-11 * tc[3]
            +5.25544750e-16 * tc[4];
        /*species 10: N2O */
        species[10] =
            +8.16581100e+03 * invT
            +5.37622700e+00
            -4.71897700e+00 * tc[0]
            -1.43685700e-03 * tc[1]
            +1.99582667e-07 * tc[2]
            -1.87546000e-11 * tc[3]
            +7.87668500e-16 * tc[4];
        /*species 11: NH */
        species[11] =
            +4.21345140e+04 * invT
            -3.95708690e+00
            -2.78369290e+00 * tc[0]
            -6.64921450e-04 * tc[1]
            +7.07967450e-08 * tc[2]
            -6.52904200e-12 * tc[3]
            +2.75222350e-16 * tc[4];
        /*species 12: N */
        species[12] =
            +5.61160400e+04 * invT
            -2.99849000e+00
            -2.45026800e+00 * tc[0]
            -5.33073000e-05 * tc[1]
            +1.24422283e-08 * tc[2]
            -1.56637667e-12 * tc[3]
            +5.12992000e-17 * tc[4];
        /*species 13: NNH */
        species[13] =
            +2.78802900e+04 * invT
            +2.51105320e+00
            -4.41534200e+00 * tc[0]
            -8.07194000e-04 * tc[1]
            +2.72149000e-08 * tc[2]
            +7.13320500e-12 * tc[3]
            -8.07395500e-16 * tc[4];
        /*species 14: N2 */
        species[14] =
            -9.23948688e+02 * invT
            -3.91931125e+00
            -2.95257637e+00 * tc[0]
            -6.98450200e-04 * tc[1]
            +8.21052672e-08 * tc[2]
            -6.55008496e-12 * tc[3]
            +2.30377602e-16 * tc[4];
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
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 2: OH */
        species[2] =
            +2.99198424e+00
            -2.40106655e-03 * tc[1]
            +4.61664033e-06 * tc[2]
            -3.87916306e-09 * tc[3]
            +1.36319502e-12 * tc[4];
        /*species 3: H2 */
        species[3] =
            +1.34430290e+00
            +7.98042480e-03 * tc[1]
            -1.94779170e-05 * tc[2]
            +2.01569670e-08 * tc[3]
            -7.37602890e-12 * tc[4];
        /*species 4: O2 */
        species[4] =
            +2.78245636e+00
            -2.99673415e-03 * tc[1]
            +9.84730200e-06 * tc[2]
            -9.68129508e-09 * tc[3]
            +3.24372836e-12 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +3.30178800e+00
            -4.74902010e-03 * tc[1]
            +2.11579530e-05 * tc[2]
            -2.42759610e-08 * tc[3]
            +9.29206700e-12 * tc[4];
        /*species 6: H2O */
        species[6] =
            +3.19863520e+00
            -2.03640170e-03 * tc[1]
            +6.52034160e-06 * tc[2]
            -5.48792690e-09 * tc[3]
            +1.77196800e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.31515149e+00
            -8.47390622e-04 * tc[1]
            +1.76404323e-05 * tc[2]
            -2.26762944e-08 * tc[3]
            +9.08950158e-12 * tc[4];
        /*species 8: NO */
        species[8] =
            +3.21859896e+00
            -4.63988124e-03 * tc[1]
            +1.10443049e-05 * tc[2]
            -9.34055507e-09 * tc[3]
            +2.80554874e-12 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +2.94403120e+00
            -1.58542900e-03 * tc[1]
            +1.66578120e-05 * tc[2]
            -2.04754260e-08 * tc[3]
            +7.83505640e-12 * tc[4];
        /*species 10: N2O */
        species[10] =
            +1.54305800e+00
            +9.49219300e-03 * tc[1]
            -9.79277500e-06 * tc[2]
            +6.26384500e-09 * tc[3]
            -1.90182600e-12 * tc[4];
        /*species 11: NH */
        species[11] =
            +2.49290840e+00
            +3.11791970e-04 * tc[1]
            -1.48904840e-06 * tc[2]
            +2.48164420e-09 * tc[3]
            -1.03569670e-12 * tc[4];
        /*species 12: N */
        species[12] =
            +1.50307100e+00
            -2.18001800e-05 * tc[1]
            +5.42052900e-08 * tc[2]
            -5.64756000e-11 * tc[3]
            +2.09990400e-14 * tc[4];
        /*species 13: NNH */
        species[13] =
            +2.50134400e+00
            +2.05358700e-03 * tc[1]
            +7.17041000e-07 * tc[2]
            +4.92134800e-10 * tc[3]
            -9.67117000e-13 * tc[4];
        /*species 14: N2 */
        species[14] =
            +2.53100528e+00
            -1.23660988e-04 * tc[1]
            -5.02999433e-07 * tc[2]
            +2.43530612e-09 * tc[3]
            -1.40881235e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +1.54363697e+00
            -2.73162486e-05 * tc[1]
            -4.19029520e-09 * tc[2]
            +4.95481845e-12 * tc[3]
            -4.79553694e-16 * tc[4];
        /*species 2: OH */
        species[2] =
            +1.83853033e+00
            +1.10741289e-03 * tc[1]
            -2.94000209e-07 * tc[2]
            +4.20698729e-11 * tc[3]
            -2.42289890e-15 * tc[4];
        /*species 3: H2 */
        species[3] =
            +1.93283050e+00
            +8.26598020e-04 * tc[1]
            -1.46400570e-07 * tc[2]
            +1.54098510e-11 * tc[3]
            -6.88796150e-16 * tc[4];
        /*species 4: O2 */
        species[4] =
            +2.66096083e+00
            +6.56365523e-04 * tc[1]
            -1.41149485e-07 * tc[2]
            +2.05797658e-11 * tc[3]
            -1.29913248e-15 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +3.17226590e+00
            +1.88120980e-03 * tc[1]
            -3.46292970e-07 * tc[2]
            +1.94685160e-11 * tc[3]
            +1.76091530e-16 * tc[4];
        /*species 6: H2O */
        species[6] =
            +1.67703890e+00
            +2.97318160e-03 * tc[1]
            -7.73768890e-07 * tc[2]
            +9.44335140e-11 * tc[3]
            -4.26899910e-15 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57977305e+00
            +4.05326003e-03 * tc[1]
            -1.29844730e-06 * tc[2]
            +1.98211400e-10 * tc[3]
            -1.13968792e-14 * tc[4];
        /*species 8: NO */
        species[8] =
            +2.26071234e+00
            +1.19101135e-03 * tc[1]
            -4.29122646e-07 * tc[2]
            +6.94481463e-11 * tc[3]
            -4.03295681e-15 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +3.88475400e+00
            +2.17239550e-03 * tc[1]
            -8.28069090e-07 * tc[2]
            +1.57475100e-10 * tc[3]
            -1.05108950e-14 * tc[4];
        /*species 10: N2O */
        species[10] =
            +3.71897700e+00
            +2.87371400e-03 * tc[1]
            -1.19749600e-06 * tc[2]
            +2.25055200e-10 * tc[3]
            -1.57533700e-14 * tc[4];
        /*species 11: NH */
        species[11] =
            +1.78369290e+00
            +1.32984290e-03 * tc[1]
            -4.24780470e-07 * tc[2]
            +7.83485040e-11 * tc[3]
            -5.50444700e-15 * tc[4];
        /*species 12: N */
        species[12] =
            +1.45026800e+00
            +1.06614600e-04 * tc[1]
            -7.46533700e-08 * tc[2]
            +1.87965200e-11 * tc[3]
            -1.02598400e-15 * tc[4];
        /*species 13: NNH */
        species[13] =
            +3.41534200e+00
            +1.61438800e-03 * tc[1]
            -1.63289400e-07 * tc[2]
            -8.55984600e-11 * tc[3]
            +1.61479100e-14 * tc[4];
        /*species 14: N2 */
        species[14] =
            +1.95257637e+00
            +1.39690040e-03 * tc[1]
            -4.92631603e-07 * tc[2]
            +7.86010195e-11 * tc[3]
            -4.60755204e-15 * tc[4];
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
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.99198424e+00
            -2.40106655e-03 * tc[1]
            +4.61664033e-06 * tc[2]
            -3.87916306e-09 * tc[3]
            +1.36319502e-12 * tc[4];
        /*species 3: H2 */
        species[3] =
            +2.34430290e+00
            +7.98042480e-03 * tc[1]
            -1.94779170e-05 * tc[2]
            +2.01569670e-08 * tc[3]
            -7.37602890e-12 * tc[4];
        /*species 4: O2 */
        species[4] =
            +3.78245636e+00
            -2.99673415e-03 * tc[1]
            +9.84730200e-06 * tc[2]
            -9.68129508e-09 * tc[3]
            +3.24372836e-12 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +4.30178800e+00
            -4.74902010e-03 * tc[1]
            +2.11579530e-05 * tc[2]
            -2.42759610e-08 * tc[3]
            +9.29206700e-12 * tc[4];
        /*species 6: H2O */
        species[6] =
            +4.19863520e+00
            -2.03640170e-03 * tc[1]
            +6.52034160e-06 * tc[2]
            -5.48792690e-09 * tc[3]
            +1.77196800e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.31515149e+00
            -8.47390622e-04 * tc[1]
            +1.76404323e-05 * tc[2]
            -2.26762944e-08 * tc[3]
            +9.08950158e-12 * tc[4];
        /*species 8: NO */
        species[8] =
            +4.21859896e+00
            -4.63988124e-03 * tc[1]
            +1.10443049e-05 * tc[2]
            -9.34055507e-09 * tc[3]
            +2.80554874e-12 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +3.94403120e+00
            -1.58542900e-03 * tc[1]
            +1.66578120e-05 * tc[2]
            -2.04754260e-08 * tc[3]
            +7.83505640e-12 * tc[4];
        /*species 10: N2O */
        species[10] =
            +2.54305800e+00
            +9.49219300e-03 * tc[1]
            -9.79277500e-06 * tc[2]
            +6.26384500e-09 * tc[3]
            -1.90182600e-12 * tc[4];
        /*species 11: NH */
        species[11] =
            +3.49290840e+00
            +3.11791970e-04 * tc[1]
            -1.48904840e-06 * tc[2]
            +2.48164420e-09 * tc[3]
            -1.03569670e-12 * tc[4];
        /*species 12: N */
        species[12] =
            +2.50307100e+00
            -2.18001800e-05 * tc[1]
            +5.42052900e-08 * tc[2]
            -5.64756000e-11 * tc[3]
            +2.09990400e-14 * tc[4];
        /*species 13: NNH */
        species[13] =
            +3.50134400e+00
            +2.05358700e-03 * tc[1]
            +7.17041000e-07 * tc[2]
            +4.92134800e-10 * tc[3]
            -9.67117000e-13 * tc[4];
        /*species 14: N2 */
        species[14] =
            +3.53100528e+00
            -1.23660988e-04 * tc[1]
            -5.02999433e-07 * tc[2]
            +2.43530612e-09 * tc[3]
            -1.40881235e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.54363697e+00
            -2.73162486e-05 * tc[1]
            -4.19029520e-09 * tc[2]
            +4.95481845e-12 * tc[3]
            -4.79553694e-16 * tc[4];
        /*species 2: OH */
        species[2] =
            +2.83853033e+00
            +1.10741289e-03 * tc[1]
            -2.94000209e-07 * tc[2]
            +4.20698729e-11 * tc[3]
            -2.42289890e-15 * tc[4];
        /*species 3: H2 */
        species[3] =
            +2.93283050e+00
            +8.26598020e-04 * tc[1]
            -1.46400570e-07 * tc[2]
            +1.54098510e-11 * tc[3]
            -6.88796150e-16 * tc[4];
        /*species 4: O2 */
        species[4] =
            +3.66096083e+00
            +6.56365523e-04 * tc[1]
            -1.41149485e-07 * tc[2]
            +2.05797658e-11 * tc[3]
            -1.29913248e-15 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +4.17226590e+00
            +1.88120980e-03 * tc[1]
            -3.46292970e-07 * tc[2]
            +1.94685160e-11 * tc[3]
            +1.76091530e-16 * tc[4];
        /*species 6: H2O */
        species[6] =
            +2.67703890e+00
            +2.97318160e-03 * tc[1]
            -7.73768890e-07 * tc[2]
            +9.44335140e-11 * tc[3]
            -4.26899910e-15 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57977305e+00
            +4.05326003e-03 * tc[1]
            -1.29844730e-06 * tc[2]
            +1.98211400e-10 * tc[3]
            -1.13968792e-14 * tc[4];
        /*species 8: NO */
        species[8] =
            +3.26071234e+00
            +1.19101135e-03 * tc[1]
            -4.29122646e-07 * tc[2]
            +6.94481463e-11 * tc[3]
            -4.03295681e-15 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +4.88475400e+00
            +2.17239550e-03 * tc[1]
            -8.28069090e-07 * tc[2]
            +1.57475100e-10 * tc[3]
            -1.05108950e-14 * tc[4];
        /*species 10: N2O */
        species[10] =
            +4.71897700e+00
            +2.87371400e-03 * tc[1]
            -1.19749600e-06 * tc[2]
            +2.25055200e-10 * tc[3]
            -1.57533700e-14 * tc[4];
        /*species 11: NH */
        species[11] =
            +2.78369290e+00
            +1.32984290e-03 * tc[1]
            -4.24780470e-07 * tc[2]
            +7.83485040e-11 * tc[3]
            -5.50444700e-15 * tc[4];
        /*species 12: N */
        species[12] =
            +2.45026800e+00
            +1.06614600e-04 * tc[1]
            -7.46533700e-08 * tc[2]
            +1.87965200e-11 * tc[3]
            -1.02598400e-15 * tc[4];
        /*species 13: NNH */
        species[13] =
            +4.41534200e+00
            +1.61438800e-03 * tc[1]
            -1.63289400e-07 * tc[2]
            -8.55984600e-11 * tc[3]
            +1.61479100e-14 * tc[4];
        /*species 14: N2 */
        species[14] =
            +2.95257637e+00
            +1.39690040e-03 * tc[1]
            -4.92631603e-07 * tc[2]
            +7.86010195e-11 * tc[3]
            -4.60755204e-15 * tc[4];
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
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 * invT;
        /*species 1: O */
        species[1] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 2: OH */
        species[2] =
            +2.99198424e+00
            -1.20053327e-03 * tc[1]
            +1.53888011e-06 * tc[2]
            -9.69790765e-10 * tc[3]
            +2.72639004e-13 * tc[4]
            +3.37165248e+03 * invT;
        /*species 3: H2 */
        species[3] =
            +1.34430290e+00
            +3.99021240e-03 * tc[1]
            -6.49263900e-06 * tc[2]
            +5.03924175e-09 * tc[3]
            -1.47520578e-12 * tc[4]
            -9.17924130e+02 * invT;
        /*species 4: O2 */
        species[4] =
            +2.78245636e+00
            -1.49836707e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745672e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 5: HO2 */
        species[5] =
            +3.30178800e+00
            -2.37451005e-03 * tc[1]
            +7.05265100e-06 * tc[2]
            -6.06899025e-09 * tc[3]
            +1.85841340e-12 * tc[4]
            +2.63190983e+02 * invT;
        /*species 6: H2O */
        species[6] =
            +3.19863520e+00
            -1.01820085e-03 * tc[1]
            +2.17344720e-06 * tc[2]
            -1.37198172e-09 * tc[3]
            +3.54393600e-13 * tc[4]
            -3.02937260e+04 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.31515149e+00
            -4.23695311e-04 * tc[1]
            +5.88014410e-06 * tc[2]
            -5.66907360e-09 * tc[3]
            +1.81790032e-12 * tc[4]
            -1.76843601e+04 * invT;
        /*species 8: NO */
        species[8] =
            +3.21859896e+00
            -2.31994062e-03 * tc[1]
            +3.68143497e-06 * tc[2]
            -2.33513877e-09 * tc[3]
            +5.61109748e-13 * tc[4]
            +9.81823786e+03 * invT;
        /*species 9: NO2 */
        species[9] =
            +2.94403120e+00
            -7.92714500e-04 * tc[1]
            +5.55260400e-06 * tc[2]
            -5.11885650e-09 * tc[3]
            +1.56701128e-12 * tc[4]
            +2.87409757e+03 * invT;
        /*species 10: N2O */
        species[10] =
            +1.54305800e+00
            +4.74609650e-03 * tc[1]
            -3.26425833e-06 * tc[2]
            +1.56596125e-09 * tc[3]
            -3.80365200e-13 * tc[4]
            +8.76510000e+03 * invT;
        /*species 11: NH */
        species[11] =
            +2.49290840e+00
            +1.55895985e-04 * tc[1]
            -4.96349467e-07 * tc[2]
            +6.20411050e-10 * tc[3]
            -2.07139340e-13 * tc[4]
            +4.18942940e+04 * invT;
        /*species 12: N */
        species[12] =
            +1.50307100e+00
            -1.09000900e-05 * tc[1]
            +1.80684300e-08 * tc[2]
            -1.41189000e-11 * tc[3]
            +4.19980800e-15 * tc[4]
            +5.60989000e+04 * invT;
        /*species 13: NNH */
        species[13] =
            +2.50134400e+00
            +1.02679350e-03 * tc[1]
            +2.39013667e-07 * tc[2]
            +1.23033700e-10 * tc[3]
            -1.93423400e-13 * tc[4]
            +2.83334700e+04 * invT;
        /*species 14: N2 */
        species[14] =
            +2.53100528e+00
            -6.18304940e-05 * tc[1]
            -1.67666478e-07 * tc[2]
            +6.08826530e-10 * tc[3]
            -2.81762470e-13 * tc[4]
            -1.04697628e+03 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 * invT;
        /*species 1: O */
        species[1] =
            +1.54363697e+00
            -1.36581243e-05 * tc[1]
            -1.39676507e-09 * tc[2]
            +1.23870461e-12 * tc[3]
            -9.59107388e-17 * tc[4]
            +2.92260120e+04 * invT;
        /*species 2: OH */
        species[2] =
            +1.83853033e+00
            +5.53706445e-04 * tc[1]
            -9.80000697e-08 * tc[2]
            +1.05174682e-11 * tc[3]
            -4.84579780e-16 * tc[4]
            +3.70056220e+03 * invT;
        /*species 3: H2 */
        species[3] =
            +1.93283050e+00
            +4.13299010e-04 * tc[1]
            -4.88001900e-08 * tc[2]
            +3.85246275e-12 * tc[3]
            -1.37759230e-16 * tc[4]
            -8.13055820e+02 * invT;
        /*species 4: O2 */
        species[4] =
            +2.66096083e+00
            +3.28182761e-04 * tc[1]
            -4.70498283e-08 * tc[2]
            +5.14494145e-12 * tc[3]
            -2.59826496e-16 * tc[4]
            -1.21597725e+03 * invT;
        /*species 5: HO2 */
        species[5] =
            +3.17226590e+00
            +9.40604900e-04 * tc[1]
            -1.15430990e-07 * tc[2]
            +4.86712900e-12 * tc[3]
            +3.52183060e-17 * tc[4]
            +3.02010736e+01 * invT;
        /*species 6: H2O */
        species[6] =
            +1.67703890e+00
            +1.48659080e-03 * tc[1]
            -2.57922963e-07 * tc[2]
            +2.36083785e-11 * tc[3]
            -8.53799820e-16 * tc[4]
            -2.98858940e+04 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.57977305e+00
            +2.02663002e-03 * tc[1]
            -4.32815767e-07 * tc[2]
            +4.95528500e-11 * tc[3]
            -2.27937584e-15 * tc[4]
            -1.79847939e+04 * invT;
        /*species 8: NO */
        species[8] =
            +2.26071234e+00
            +5.95505675e-04 * tc[1]
            -1.43040882e-07 * tc[2]
            +1.73620366e-11 * tc[3]
            -8.06591362e-16 * tc[4]
            +9.89456954e+03 * invT;
        /*species 9: NO2 */
        species[9] =
            +3.88475400e+00
            +1.08619775e-03 * tc[1]
            -2.76023030e-07 * tc[2]
            +3.93687750e-11 * tc[3]
            -2.10217900e-15 * tc[4]
            +2.29397777e+03 * invT;
        /*species 10: N2O */
        species[10] =
            +3.71897700e+00
            +1.43685700e-03 * tc[1]
            -3.99165333e-07 * tc[2]
            +5.62638000e-11 * tc[3]
            -3.15067400e-15 * tc[4]
            +8.16581100e+03 * invT;
        /*species 11: NH */
        species[11] =
            +1.78369290e+00
            +6.64921450e-04 * tc[1]
            -1.41593490e-07 * tc[2]
            +1.95871260e-11 * tc[3]
            -1.10088940e-15 * tc[4]
            +4.21345140e+04 * invT;
        /*species 12: N */
        species[12] =
            +1.45026800e+00
            +5.33073000e-05 * tc[1]
            -2.48844567e-08 * tc[2]
            +4.69913000e-12 * tc[3]
            -2.05196800e-16 * tc[4]
            +5.61160400e+04 * invT;
        /*species 13: NNH */
        species[13] =
            +3.41534200e+00
            +8.07194000e-04 * tc[1]
            -5.44298000e-08 * tc[2]
            -2.13996150e-11 * tc[3]
            +3.22958200e-15 * tc[4]
            +2.78802900e+04 * invT;
        /*species 14: N2 */
        species[14] =
            +1.95257637e+00
            +6.98450200e-04 * tc[1]
            -1.64210534e-07 * tc[2]
            +1.96502549e-11 * tc[3]
            -9.21510408e-16 * tc[4]
            -9.23948688e+02 * invT;
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
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 * invT;
        /*species 1: O */
        species[1] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 2: OH */
        species[2] =
            +3.99198424e+00
            -1.20053327e-03 * tc[1]
            +1.53888011e-06 * tc[2]
            -9.69790765e-10 * tc[3]
            +2.72639004e-13 * tc[4]
            +3.37165248e+03 * invT;
        /*species 3: H2 */
        species[3] =
            +2.34430290e+00
            +3.99021240e-03 * tc[1]
            -6.49263900e-06 * tc[2]
            +5.03924175e-09 * tc[3]
            -1.47520578e-12 * tc[4]
            -9.17924130e+02 * invT;
        /*species 4: O2 */
        species[4] =
            +3.78245636e+00
            -1.49836707e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745672e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 5: HO2 */
        species[5] =
            +4.30178800e+00
            -2.37451005e-03 * tc[1]
            +7.05265100e-06 * tc[2]
            -6.06899025e-09 * tc[3]
            +1.85841340e-12 * tc[4]
            +2.63190983e+02 * invT;
        /*species 6: H2O */
        species[6] =
            +4.19863520e+00
            -1.01820085e-03 * tc[1]
            +2.17344720e-06 * tc[2]
            -1.37198172e-09 * tc[3]
            +3.54393600e-13 * tc[4]
            -3.02937260e+04 * invT;
        /*species 7: H2O2 */
        species[7] =
            +4.31515149e+00
            -4.23695311e-04 * tc[1]
            +5.88014410e-06 * tc[2]
            -5.66907360e-09 * tc[3]
            +1.81790032e-12 * tc[4]
            -1.76843601e+04 * invT;
        /*species 8: NO */
        species[8] =
            +4.21859896e+00
            -2.31994062e-03 * tc[1]
            +3.68143497e-06 * tc[2]
            -2.33513877e-09 * tc[3]
            +5.61109748e-13 * tc[4]
            +9.81823786e+03 * invT;
        /*species 9: NO2 */
        species[9] =
            +3.94403120e+00
            -7.92714500e-04 * tc[1]
            +5.55260400e-06 * tc[2]
            -5.11885650e-09 * tc[3]
            +1.56701128e-12 * tc[4]
            +2.87409757e+03 * invT;
        /*species 10: N2O */
        species[10] =
            +2.54305800e+00
            +4.74609650e-03 * tc[1]
            -3.26425833e-06 * tc[2]
            +1.56596125e-09 * tc[3]
            -3.80365200e-13 * tc[4]
            +8.76510000e+03 * invT;
        /*species 11: NH */
        species[11] =
            +3.49290840e+00
            +1.55895985e-04 * tc[1]
            -4.96349467e-07 * tc[2]
            +6.20411050e-10 * tc[3]
            -2.07139340e-13 * tc[4]
            +4.18942940e+04 * invT;
        /*species 12: N */
        species[12] =
            +2.50307100e+00
            -1.09000900e-05 * tc[1]
            +1.80684300e-08 * tc[2]
            -1.41189000e-11 * tc[3]
            +4.19980800e-15 * tc[4]
            +5.60989000e+04 * invT;
        /*species 13: NNH */
        species[13] =
            +3.50134400e+00
            +1.02679350e-03 * tc[1]
            +2.39013667e-07 * tc[2]
            +1.23033700e-10 * tc[3]
            -1.93423400e-13 * tc[4]
            +2.83334700e+04 * invT;
        /*species 14: N2 */
        species[14] =
            +3.53100528e+00
            -6.18304940e-05 * tc[1]
            -1.67666478e-07 * tc[2]
            +6.08826530e-10 * tc[3]
            -2.81762470e-13 * tc[4]
            -1.04697628e+03 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 * invT;
        /*species 1: O */
        species[1] =
            +2.54363697e+00
            -1.36581243e-05 * tc[1]
            -1.39676507e-09 * tc[2]
            +1.23870461e-12 * tc[3]
            -9.59107388e-17 * tc[4]
            +2.92260120e+04 * invT;
        /*species 2: OH */
        species[2] =
            +2.83853033e+00
            +5.53706445e-04 * tc[1]
            -9.80000697e-08 * tc[2]
            +1.05174682e-11 * tc[3]
            -4.84579780e-16 * tc[4]
            +3.70056220e+03 * invT;
        /*species 3: H2 */
        species[3] =
            +2.93283050e+00
            +4.13299010e-04 * tc[1]
            -4.88001900e-08 * tc[2]
            +3.85246275e-12 * tc[3]
            -1.37759230e-16 * tc[4]
            -8.13055820e+02 * invT;
        /*species 4: O2 */
        species[4] =
            +3.66096083e+00
            +3.28182761e-04 * tc[1]
            -4.70498283e-08 * tc[2]
            +5.14494145e-12 * tc[3]
            -2.59826496e-16 * tc[4]
            -1.21597725e+03 * invT;
        /*species 5: HO2 */
        species[5] =
            +4.17226590e+00
            +9.40604900e-04 * tc[1]
            -1.15430990e-07 * tc[2]
            +4.86712900e-12 * tc[3]
            +3.52183060e-17 * tc[4]
            +3.02010736e+01 * invT;
        /*species 6: H2O */
        species[6] =
            +2.67703890e+00
            +1.48659080e-03 * tc[1]
            -2.57922963e-07 * tc[2]
            +2.36083785e-11 * tc[3]
            -8.53799820e-16 * tc[4]
            -2.98858940e+04 * invT;
        /*species 7: H2O2 */
        species[7] =
            +4.57977305e+00
            +2.02663002e-03 * tc[1]
            -4.32815767e-07 * tc[2]
            +4.95528500e-11 * tc[3]
            -2.27937584e-15 * tc[4]
            -1.79847939e+04 * invT;
        /*species 8: NO */
        species[8] =
            +3.26071234e+00
            +5.95505675e-04 * tc[1]
            -1.43040882e-07 * tc[2]
            +1.73620366e-11 * tc[3]
            -8.06591362e-16 * tc[4]
            +9.89456954e+03 * invT;
        /*species 9: NO2 */
        species[9] =
            +4.88475400e+00
            +1.08619775e-03 * tc[1]
            -2.76023030e-07 * tc[2]
            +3.93687750e-11 * tc[3]
            -2.10217900e-15 * tc[4]
            +2.29397777e+03 * invT;
        /*species 10: N2O */
        species[10] =
            +4.71897700e+00
            +1.43685700e-03 * tc[1]
            -3.99165333e-07 * tc[2]
            +5.62638000e-11 * tc[3]
            -3.15067400e-15 * tc[4]
            +8.16581100e+03 * invT;
        /*species 11: NH */
        species[11] =
            +2.78369290e+00
            +6.64921450e-04 * tc[1]
            -1.41593490e-07 * tc[2]
            +1.95871260e-11 * tc[3]
            -1.10088940e-15 * tc[4]
            +4.21345140e+04 * invT;
        /*species 12: N */
        species[12] =
            +2.45026800e+00
            +5.33073000e-05 * tc[1]
            -2.48844567e-08 * tc[2]
            +4.69913000e-12 * tc[3]
            -2.05196800e-16 * tc[4]
            +5.61160400e+04 * invT;
        /*species 13: NNH */
        species[13] =
            +4.41534200e+00
            +8.07194000e-04 * tc[1]
            -5.44298000e-08 * tc[2]
            -2.13996150e-11 * tc[3]
            +3.22958200e-15 * tc[4]
            +2.78802900e+04 * invT;
        /*species 14: N2 */
        species[14] =
            +2.95257637e+00
            +6.98450200e-04 * tc[1]
            -1.64210534e-07 * tc[2]
            +1.96502549e-11 * tc[3]
            -9.21510408e-16 * tc[4]
            -9.23948688e+02 * invT;
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
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.46682850e-01 ;
        /*species 1: O */
        species[1] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 2: OH */
        species[2] =
            +3.99198424e+00 * tc[0]
            -2.40106655e-03 * tc[1]
            +2.30832017e-06 * tc[2]
            -1.29305435e-09 * tc[3]
            +3.40798755e-13 * tc[4]
            -1.03814059e-01 ;
        /*species 3: H2 */
        species[3] =
            +2.34430290e+00 * tc[0]
            +7.98042480e-03 * tc[1]
            -9.73895850e-06 * tc[2]
            +6.71898900e-09 * tc[3]
            -1.84400722e-12 * tc[4]
            +6.83002180e-01 ;
        /*species 4: O2 */
        species[4] =
            +3.78245636e+00 * tc[0]
            -2.99673415e-03 * tc[1]
            +4.92365100e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932090e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 5: HO2 */
        species[5] =
            +4.30178800e+00 * tc[0]
            -4.74902010e-03 * tc[1]
            +1.05789765e-05 * tc[2]
            -8.09198700e-09 * tc[3]
            +2.32301675e-12 * tc[4]
            +3.71587740e+00 ;
        /*species 6: H2O */
        species[6] =
            +4.19863520e+00 * tc[0]
            -2.03640170e-03 * tc[1]
            +3.26017080e-06 * tc[2]
            -1.82930897e-09 * tc[3]
            +4.42992000e-13 * tc[4]
            -8.49009010e-01 ;
        /*species 7: H2O2 */
        species[7] =
            +4.31515149e+00 * tc[0]
            -8.47390622e-04 * tc[1]
            +8.82021615e-06 * tc[2]
            -7.55876480e-09 * tc[3]
            +2.27237539e-12 * tc[4]
            +3.27373216e+00 ;
        /*species 8: NO */
        species[8] =
            +4.21859896e+00 * tc[0]
            -4.63988124e-03 * tc[1]
            +5.52215245e-06 * tc[2]
            -3.11351836e-09 * tc[3]
            +7.01387185e-13 * tc[4]
            +2.28060952e+00 ;
        /*species 9: NO2 */
        species[9] =
            +3.94403120e+00 * tc[0]
            -1.58542900e-03 * tc[1]
            +8.32890600e-06 * tc[2]
            -6.82514200e-09 * tc[3]
            +1.95876410e-12 * tc[4]
            +6.31199190e+00 ;
        /*species 10: N2O */
        species[10] =
            +2.54305800e+00 * tc[0]
            +9.49219300e-03 * tc[1]
            -4.89638750e-06 * tc[2]
            +2.08794833e-09 * tc[3]
            -4.75456500e-13 * tc[4]
            +9.51122200e+00 ;
        /*species 11: NH */
        species[11] =
            +3.49290840e+00 * tc[0]
            +3.11791970e-04 * tc[1]
            -7.44524200e-07 * tc[2]
            +8.27214733e-10 * tc[3]
            -2.58924175e-13 * tc[4]
            +1.84832770e+00 ;
        /*species 12: N */
        species[12] =
            +2.50307100e+00 * tc[0]
            -2.18001800e-05 * tc[1]
            +2.71026450e-08 * tc[2]
            -1.88252000e-11 * tc[3]
            +5.24976000e-15 * tc[4]
            +4.16756600e+00 ;
        /*species 13: NNH */
        species[13] =
            +3.50134400e+00 * tc[0]
            +2.05358700e-03 * tc[1]
            +3.58520500e-07 * tc[2]
            +1.64044933e-10 * tc[3]
            -2.41779250e-13 * tc[4]
            +6.39183700e+00 ;
        /*species 14: N2 */
        species[14] =
            +3.53100528e+00 * tc[0]
            -1.23660988e-04 * tc[1]
            -2.51499717e-07 * tc[2]
            +8.11768707e-10 * tc[3]
            -3.52203088e-13 * tc[4]
            +2.96747038e+00 ;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.46682850e-01 ;
        /*species 1: O */
        species[1] =
            +2.54363697e+00 * tc[0]
            -2.73162486e-05 * tc[1]
            -2.09514760e-09 * tc[2]
            +1.65160615e-12 * tc[3]
            -1.19888423e-16 * tc[4]
            +4.92229457e+00 ;
        /*species 2: OH */
        species[2] =
            +2.83853033e+00 * tc[0]
            +1.10741289e-03 * tc[1]
            -1.47000104e-07 * tc[2]
            +1.40232910e-11 * tc[3]
            -6.05724725e-16 * tc[4]
            +5.84513094e+00 ;
        /*species 3: H2 */
        species[3] =
            +2.93283050e+00 * tc[0]
            +8.26598020e-04 * tc[1]
            -7.32002850e-08 * tc[2]
            +5.13661700e-12 * tc[3]
            -1.72199037e-16 * tc[4]
            -1.02431640e+00 ;
        /*species 4: O2 */
        species[4] =
            +3.66096083e+00 * tc[0]
            +6.56365523e-04 * tc[1]
            -7.05747425e-08 * tc[2]
            +6.85992193e-12 * tc[3]
            -3.24783120e-16 * tc[4]
            +3.41536184e+00 ;
        /*species 5: HO2 */
        species[5] =
            +4.17226590e+00 * tc[0]
            +1.88120980e-03 * tc[1]
            -1.73146485e-07 * tc[2]
            +6.48950533e-12 * tc[3]
            +4.40228825e-17 * tc[4]
            +2.95697380e+00 ;
        /*species 6: H2O */
        species[6] =
            +2.67703890e+00 * tc[0]
            +2.97318160e-03 * tc[1]
            -3.86884445e-07 * tc[2]
            +3.14778380e-11 * tc[3]
            -1.06724977e-15 * tc[4]
            +6.88255000e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.57977305e+00 * tc[0]
            +4.05326003e-03 * tc[1]
            -6.49223650e-07 * tc[2]
            +6.60704667e-11 * tc[3]
            -2.84921980e-15 * tc[4]
            +6.64969660e-01 ;
        /*species 8: NO */
        species[8] =
            +3.26071234e+00 * tc[0]
            +1.19101135e-03 * tc[1]
            -2.14561323e-07 * tc[2]
            +2.31493821e-11 * tc[3]
            -1.00823920e-15 * tc[4]
            +6.36900469e+00 ;
        /*species 9: NO2 */
        species[9] =
            +4.88475400e+00 * tc[0]
            +2.17239550e-03 * tc[1]
            -4.14034545e-07 * tc[2]
            +5.24917000e-11 * tc[3]
            -2.62772375e-15 * tc[4]
            -1.17416951e-01 ;
        /*species 10: N2O */
        species[10] =
            +4.71897700e+00 * tc[0]
            +2.87371400e-03 * tc[1]
            -5.98748000e-07 * tc[2]
            +7.50184000e-11 * tc[3]
            -3.93834250e-15 * tc[4]
            -1.65725000e+00 ;
        /*species 11: NH */
        species[11] =
            +2.78369290e+00 * tc[0]
            +1.32984290e-03 * tc[1]
            -2.12390235e-07 * tc[2]
            +2.61161680e-11 * tc[3]
            -1.37611175e-15 * tc[4]
            +5.74077980e+00 ;
        /*species 12: N */
        species[12] =
            +2.45026800e+00 * tc[0]
            +1.06614600e-04 * tc[1]
            -3.73266850e-08 * tc[2]
            +6.26550667e-12 * tc[3]
            -2.56496000e-16 * tc[4]
            +4.44875800e+00 ;
        /*species 13: NNH */
        species[13] =
            +4.41534200e+00 * tc[0]
            +1.61438800e-03 * tc[1]
            -8.16447000e-08 * tc[2]
            -2.85328200e-11 * tc[3]
            +4.03697750e-15 * tc[4]
            +9.04288800e-01 ;
        /*species 14: N2 */
        species[14] =
            +2.95257637e+00 * tc[0]
            +1.39690040e-03 * tc[1]
            -2.46315802e-07 * tc[2]
            +2.62003398e-11 * tc[3]
            -1.15188801e-15 * tc[4]
            +5.87188762e+00 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * restrict wt)
{
    wt[0] = 1.007970; /*H */
    wt[1] = 15.999400; /*O */
    wt[2] = 17.007370; /*OH */
    wt[3] = 2.015940; /*H2 */
    wt[4] = 31.998800; /*O2 */
    wt[5] = 33.006770; /*HO2 */
    wt[6] = 18.015340; /*H2O */
    wt[7] = 34.014740; /*H2O2 */
    wt[8] = 30.006100; /*NO */
    wt[9] = 46.005500; /*NO2 */
    wt[10] = 44.012800; /*N2O */
    wt[11] = 15.014670; /*NH */
    wt[12] = 14.006700; /*N */
    wt[13] = 29.021370; /*NNH */
    wt[14] = 28.013400; /*N2 */

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
  *LENIMC =           62;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC) {
  *LENRMC =         4800;}
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
  *KK =           15;}
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
  WT[           0] =   0.1007969975471497E+01;
  WT[           1] =   0.1599940013885498E+02;
  WT[           2] =   0.1700737011432648E+02;
  WT[           3] =   0.2015939950942993E+01;
  WT[           4] =   0.3199880027770996E+02;
  WT[           5] =   0.3300677025318146E+02;
  WT[           6] =   0.1801534008979797E+02;
  WT[           7] =   0.3401474022865295E+02;
  WT[           8] =   0.3000609970092773E+02;
  WT[           9] =   0.4600549983978271E+02;
  WT[          10] =   0.4401279926300049E+02;
  WT[          11] =   0.1501466953754425E+02;
  WT[          12] =   0.1400669956207275E+02;
  WT[          13] =   0.2902136909961700E+02;
  WT[          14] =   0.2801339912414551E+02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
void egtransetEPS(double* EPS) {
  EPS[           0] =   0.1450000000000000E+03;
  EPS[           1] =   0.8000000000000000E+02;
  EPS[           2] =   0.8000000000000000E+02;
  EPS[           3] =   0.3800000000000000E+02;
  EPS[           4] =   0.1074000000000000E+03;
  EPS[           5] =   0.1074000000000000E+03;
  EPS[           6] =   0.5724000000000000E+03;
  EPS[           7] =   0.1074000000000000E+03;
  EPS[           8] =   0.9753000000000000E+02;
  EPS[           9] =   0.2000000000000000E+03;
  EPS[          10] =   0.2324000000000000E+03;
  EPS[          11] =   0.8000000000000000E+02;
  EPS[          12] =   0.7140000000000001E+02;
  EPS[          13] =   0.7140000000000001E+02;
  EPS[          14] =   0.9753000000000000E+02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void egtransetSIG(double* SIG) {
  SIG[           0] =   0.2050000000000000E+01;
  SIG[           1] =   0.2750000000000000E+01;
  SIG[           2] =   0.2750000000000000E+01;
  SIG[           3] =   0.2920000000000000E+01;
  SIG[           4] =   0.3458000000000000E+01;
  SIG[           5] =   0.3458000000000000E+01;
  SIG[           6] =   0.2605000000000000E+01;
  SIG[           7] =   0.3458000000000000E+01;
  SIG[           8] =   0.3621000000000000E+01;
  SIG[           9] =   0.3500000000000000E+01;
  SIG[          10] =   0.3828000000000000E+01;
  SIG[          11] =   0.2650000000000000E+01;
  SIG[          12] =   0.3298000000000000E+01;
  SIG[          13] =   0.3798000000000000E+01;
  SIG[          14] =   0.3621000000000000E+01;
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
  DIP[           6] =   0.1844000000000000E+01;
  DIP[           7] =   0.0000000000000000E+00;
  DIP[           8] =   0.0000000000000000E+00;
  DIP[           9] =   0.0000000000000000E+00;
  DIP[          10] =   0.0000000000000000E+00;
  DIP[          11] =   0.0000000000000000E+00;
  DIP[          12] =   0.0000000000000000E+00;
  DIP[          13] =   0.0000000000000000E+00;
  DIP[          14] =   0.0000000000000000E+00;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPOL EGTRANSETPOL
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPOL egtransetpol
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPOL egtransetpol_
#endif
void egtransetPOL(double* POL) {
  POL[           0] =   0.0000000000000000E+00;
  POL[           1] =   0.0000000000000000E+00;
  POL[           2] =   0.0000000000000000E+00;
  POL[           3] =   0.7900000000000000E+00;
  POL[           4] =   0.1600000000000000E+01;
  POL[           5] =   0.0000000000000000E+00;
  POL[           6] =   0.0000000000000000E+00;
  POL[           7] =   0.0000000000000000E+00;
  POL[           8] =   0.1760000000000000E+01;
  POL[           9] =   0.0000000000000000E+00;
  POL[          10] =   0.0000000000000000E+00;
  POL[          11] =   0.0000000000000000E+00;
  POL[          12] =   0.0000000000000000E+00;
  POL[          13] =   0.0000000000000000E+00;
  POL[          14] =   0.1760000000000000E+01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetZROT EGTRANSETZROT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetZROT egtransetzrot
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetZROT egtransetzrot_
#endif
void egtransetZROT(double* ZROT) {
  ZROT[           0] =   0.0000000000000000E+00;
  ZROT[           1] =   0.0000000000000000E+00;
  ZROT[           2] =   0.0000000000000000E+00;
  ZROT[           3] =   0.2800000000000000E+03;
  ZROT[           4] =   0.3800000000000000E+01;
  ZROT[           5] =   0.1000000000000000E+01;
  ZROT[           6] =   0.4000000000000000E+01;
  ZROT[           7] =   0.3800000000000000E+01;
  ZROT[           8] =   0.4000000000000000E+01;
  ZROT[           9] =   0.1000000000000000E+01;
  ZROT[          10] =   0.1000000000000000E+01;
  ZROT[          11] =   0.4000000000000000E+01;
  ZROT[          12] =   0.0000000000000000E+00;
  ZROT[          13] =   0.1000000000000000E+01;
  ZROT[          14] =   0.4000000000000000E+01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLIN EGTRANSETNLIN
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLIN egtransetnlin
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLIN egtransetnlin_
#endif
void egtransetNLIN(int* NLIN) {
  NLIN[           0] =            0;
  NLIN[           1] =            0;
  NLIN[           2] =            1;
  NLIN[           3] =            1;
  NLIN[           4] =            1;
  NLIN[           5] =            2;
  NLIN[           6] =            2;
  NLIN[           7] =            2;
  NLIN[           8] =            1;
  NLIN[           9] =            2;
  NLIN[          10] =            1;
  NLIN[          11] =            1;
  NLIN[          12] =            0;
  NLIN[          13] =            2;
  NLIN[          14] =            1;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
  COFLAM[           0] =   0.2275028008756871E-01;
  COFLAM[           1] =   0.3261057815873756E+01;
  COFLAM[           2] =  -0.3401775687039983E+00;
  COFLAM[           3] =   0.1473938515985838E-01;
  COFLAM[           4] =   0.2241322383701714E+01;
  COFLAM[           5] =   0.1680892722599026E+01;
  COFLAM[           6] =  -0.1372174737383017E+00;
  COFLAM[           7] =   0.6048682074622076E-02;
  COFLAM[           8] =   0.1512304828773671E+02;
  COFLAM[           9] =  -0.3689219490819916E+01;
  COFLAM[          10] =   0.6016560270677763E+00;
  COFLAM[          11] =  -0.2671941052088746E-01;
  COFLAM[          12] =   0.1216294690770399E+02;
  COFLAM[          13] =  -0.1793276805201434E+01;
  COFLAM[          14] =   0.3145187605615077E+00;
  COFLAM[          15] =  -0.1247107803117781E-01;
  COFLAM[          16] =  -0.2574697546797322E+01;
  COFLAM[          17] =   0.3192053452600490E+01;
  COFLAM[          18] =  -0.3180048007278750E+00;
  COFLAM[          19] =   0.1394841585373175E-01;
  COFLAM[          20] =   0.2445641266014231E+01;
  COFLAM[          21] =   0.7562095410461581E+00;
  COFLAM[          22] =   0.6924571125398163E-01;
  COFLAM[          23] =  -0.5505326490175595E-02;
  COFLAM[          24] =   0.2078007475928601E+02;
  COFLAM[          25] =  -0.7829322376892079E+01;
  COFLAM[          26] =   0.1363321975831565E+01;
  COFLAM[          27] =  -0.6791541865700287E-01;
  COFLAM[          28] =   0.6967028474928928E+00;
  COFLAM[          29] =   0.1408219306022218E+01;
  COFLAM[          30] =   0.4097067816780268E-02;
  COFLAM[          31] =  -0.3613998188093195E-02;
  COFLAM[          32] =   0.6945977840354891E+01;
  COFLAM[          33] =  -0.9638832761755233E+00;
  COFLAM[          34] =   0.2821741717745234E+00;
  COFLAM[          35] =  -0.1493044672153078E-01;
  COFLAM[          36] =  -0.1662184858337771E+02;
  COFLAM[          37] =   0.8484044413292917E+01;
  COFLAM[          38] =  -0.9781359188964711E+00;
  COFLAM[          39] =   0.4140341324555693E-01;
  COFLAM[          40] =  -0.1212887626987692E+02;
  COFLAM[          41] =   0.6295435044043567E+01;
  COFLAM[          42] =  -0.6377803898331278E+00;
  COFLAM[          43] =   0.2394195295432890E-01;
  COFLAM[          44] =   0.1041449118465093E+02;
  COFLAM[          45] =  -0.1636552719768927E+01;
  COFLAM[          46] =   0.3162470266930527E+00;
  COFLAM[          47] =  -0.1346365879961693E-01;
  COFLAM[          48] =   0.2571501390604503E+01;
  COFLAM[          49] =   0.1435550225797724E+01;
  COFLAM[          50] =  -0.1043476296606175E+00;
  COFLAM[          51] =   0.4580649134638923E-02;
  COFLAM[          52] =   0.4921299057063938E+01;
  COFLAM[          53] =  -0.4374718231394432E+00;
  COFLAM[          54] =   0.2592263777073271E+00;
  COFLAM[          55] =  -0.1567139518196203E-01;
  COFLAM[          56] =   0.1069769027684588E+02;
  COFLAM[          57] =  -0.2529619934445385E+01;
  COFLAM[          58] =   0.4974738094590221E+00;
  COFLAM[          59] =  -0.2466468043187436E-01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void egtransetCOFETA(double* COFETA) {
  COFETA[           0] =  -0.1952716500914950E+02;
  COFETA[           1] =   0.3261057815873812E+01;
  COFETA[           2] =  -0.3401775687040061E+00;
  COFETA[           3] =   0.1473938515985874E-01;
  COFETA[           4] =  -0.1454398005828718E+02;
  COFETA[           5] =   0.1680892722599193E+01;
  COFETA[           6] =  -0.1372174737383249E+00;
  COFETA[           7] =   0.6048682074623137E-02;
  COFETA[           8] =  -0.1451343228027313E+02;
  COFETA[           9] =   0.1680892722599180E+01;
  COFETA[          10] =  -0.1372174737383228E+00;
  COFETA[          11] =   0.6048682074623037E-02;
  COFETA[          12] =  -0.1350232180202826E+02;
  COFETA[          13] =   0.8541080576610959E+00;
  COFETA[          14] =  -0.2829855499940444E-01;
  COFETA[          15] =   0.1271302744097607E-02;
  COFETA[          16] =  -0.1652175091679111E+02;
  COFETA[          17] =   0.2394046741797939E+01;
  COFETA[          18] =  -0.2301637580705922E+00;
  COFETA[          19] =   0.1008479172141301E-01;
  COFETA[          20] =  -0.1650624377237309E+02;
  COFETA[          21] =   0.2394046741797971E+01;
  COFETA[          22] =  -0.2301637580705965E+00;
  COFETA[          23] =   0.1008479172141321E-01;
  COFETA[          24] =  -0.1307160161390482E+02;
  COFETA[          25] =  -0.2592861555748340E+00;
  COFETA[          26] =   0.2565297612543304E+00;
  COFETA[          27] =  -0.1609270043920451E-01;
  COFETA[          28] =  -0.1649120313877725E+02;
  COFETA[          29] =   0.2394046741798005E+01;
  COFETA[          30] =  -0.2301637580706010E+00;
  COFETA[          31] =   0.1008479172141341E-01;
  COFETA[          32] =  -0.1595814436337888E+02;
  COFETA[          33] =   0.2131939570967141E+01;
  COFETA[          34] =  -0.1961911683803149E+00;
  COFETA[          35] =   0.8617691675965194E-02;
  COFETA[          36] =  -0.2156507583481173E+02;
  COFETA[          37] =   0.4305227753578706E+01;
  COFETA[          38] =  -0.4687113619224428E+00;
  COFETA[          39] =   0.2002099270261577E-01;
  COFETA[          40] =  -0.2301694015729145E+02;
  COFETA[          41] =   0.4733575530530014E+01;
  COFETA[          42] =  -0.5183508824318821E+00;
  COFETA[          43] =   0.2193212281602507E-01;
  COFETA[          44] =  -0.1450165928384816E+02;
  COFETA[          45] =   0.1680892722599088E+01;
  COFETA[          46] =  -0.1372174737383104E+00;
  COFETA[          47] =   0.6048682074622478E-02;
  COFETA[          48] =  -0.1434681652629993E+02;
  COFETA[          49] =   0.1435550225797739E+01;
  COFETA[          50] =  -0.1043476296606198E+00;
  COFETA[          51] =   0.4580649134639032E-02;
  COFETA[          52] =  -0.1426488496838577E+02;
  COFETA[          53] =   0.1435550225797706E+01;
  COFETA[          54] =  -0.1043476296606149E+00;
  COFETA[          55] =   0.4580649134638796E-02;
  COFETA[          56] =  -0.1599250323762854E+02;
  COFETA[          57] =   0.2131939570967079E+01;
  COFETA[          58] =  -0.1961911683803056E+00;
  COFETA[          59] =   0.8617691675964734E-02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void egtransetCOFD(double* COFD) {
  COFD[           0] =  -0.1453104126392958E+02;
  COFD[           1] =   0.4090832270068148E+01;
  COFD[           2] =  -0.3125064363694409E+00;
  COFD[           3] =   0.1336815462999992E-01;
  COFD[           4] =  -0.1446347541294260E+02;
  COFD[           5] =   0.3893324734859439E+01;
  COFD[           6] =  -0.2926671223589983E+00;
  COFD[           7] =   0.1274254538066624E-01;
  COFD[           8] =  -0.1447793881138075E+02;
  COFD[           9] =   0.3898507010768475E+01;
  COFD[          10] =  -0.2933558378767546E+00;
  COFD[          11] =   0.1277297948985323E-01;
  COFD[          12] =  -0.1127823399275531E+02;
  COFD[          13] =   0.2707816433593029E+01;
  COFD[          14] =  -0.1386141258544173E+00;
  COFD[          15] =   0.6070036436918267E-02;
  COFD[          16] =  -0.1642636080356586E+02;
  COFD[          17] =   0.4526058244629064E+01;
  COFD[          18] =  -0.3728782124089335E+00;
  COFD[          19] =   0.1613117546672936E-01;
  COFD[          20] =  -0.1643428330282377E+02;
  COFD[          21] =   0.4529060545715668E+01;
  COFD[          22] =  -0.3732727663838795E+00;
  COFD[          23] =   0.1614842241690972E-01;
  COFD[          24] =  -0.1734880949753627E+02;
  COFD[          25] =   0.4593172772581599E+01;
  COFD[          26] =  -0.3376039007204080E+00;
  COFD[          27] =   0.1287220770047757E-01;
  COFD[          28] =  -0.1644185676789823E+02;
  COFD[          29] =   0.4531933357015490E+01;
  COFD[          30] =  -0.3736503142227166E+00;
  COFD[          31] =   0.1616492637429652E-01;
  COFD[          32] =  -0.1655778410328815E+02;
  COFD[          33] =   0.4581623131684943E+01;
  COFD[          34] =  -0.3829686886254905E+00;
  COFD[          35] =   0.1668786826439446E-01;
  COFD[          36] =  -0.1765762436087865E+02;
  COFD[          37] =   0.4866441430895612E+01;
  COFD[          38] =  -0.4030149143618096E+00;
  COFD[          39] =   0.1685813093438678E-01;
  COFD[          40] =  -0.1786291016737098E+02;
  COFD[          41] =   0.4848948909820512E+01;
  COFD[          42] =  -0.3945292942419953E+00;
  COFD[          43] =   0.1623802827748926E-01;
  COFD[          44] =  -0.1434851661908174E+02;
  COFD[          45] =   0.3864395688187810E+01;
  COFD[          46] =  -0.2888220557172488E+00;
  COFD[          47] =   0.1257261544647088E-01;
  COFD[          48] =  -0.1473052091803629E+02;
  COFD[          49] =   0.3938367747699902E+01;
  COFD[          50] =  -0.3007580570513714E+00;
  COFD[          51] =   0.1318674711212272E-01;
  COFD[          52] =  -0.1555514842196157E+02;
  COFD[          53] =   0.4198163984729581E+01;
  COFD[          54] =  -0.3357679563943392E+00;
  COFD[          55] =   0.1475425817654628E-01;
  COFD[          56] =  -0.1653528873346560E+02;
  COFD[          57] =   0.4572908957761951E+01;
  COFD[          58] =  -0.3818057701602784E+00;
  COFD[          59] =   0.1663628005999293E-01;
  COFD[          60] =  -0.1446347541294260E+02;
  COFD[          61] =   0.3893324734859439E+01;
  COFD[          62] =  -0.2926671223589983E+00;
  COFD[          63] =   0.1274254538066624E-01;
  COFD[          64] =  -0.1303230661459400E+02;
  COFD[          65] =   0.2821999540521930E+01;
  COFD[          66] =  -0.1533715125660948E+00;
  COFD[          67] =   0.6707714606075767E-02;
  COFD[          68] =  -0.1304821435550623E+02;
  COFD[          69] =   0.2822365020409423E+01;
  COFD[          70] =  -0.1534211775973923E+00;
  COFD[          71] =   0.6709953971249074E-02;
  COFD[          72] =  -0.1059938809341299E+02;
  COFD[          73] =   0.2155942578374456E+01;
  COFD[          74] =  -0.6501315584613947E-01;
  COFD[          75] =   0.2795438090048244E-02;
  COFD[          76] =  -0.1443514465604553E+02;
  COFD[          77] =   0.3212703400637984E+01;
  COFD[          78] =  -0.2043100955529739E+00;
  COFD[          79] =   0.8921339110627811E-02;
  COFD[          80] =  -0.1445309446773090E+02;
  COFD[          81] =   0.3218021790448843E+01;
  COFD[          82] =  -0.2050253330481837E+00;
  COFD[          83] =   0.8953286991448825E-02;
  COFD[          84] =  -0.1876832382228059E+02;
  COFD[          85] =   0.4857884811542586E+01;
  COFD[          86] =  -0.4000862896404223E+00;
  COFD[          87] =   0.1668799943730030E-01;
  COFD[          88] =  -0.1447055666380567E+02;
  COFD[          89] =   0.3223243595288925E+01;
  COFD[          90] =  -0.2057275599882618E+00;
  COFD[          91] =   0.8984653058524278E-02;
  COFD[          92] =  -0.1421921207753527E+02;
  COFD[          93] =   0.3119581660542013E+01;
  COFD[          94] =  -0.1927352195932762E+00;
  COFD[          95] =   0.8442234647456850E-02;
  COFD[          96] =  -0.1666875714366874E+02;
  COFD[          97] =   0.4017058659430817E+01;
  COFD[          98] =  -0.3059783073756657E+00;
  COFD[          99] =   0.1320711908953701E-01;
  COFD[         100] =  -0.1720382549451406E+02;
  COFD[         101] =   0.4167339155910690E+01;
  COFD[         102] =  -0.3234483507589772E+00;
  COFD[         103] =   0.1387892452356678E-01;
  COFD[         104] =  -0.1298157331946940E+02;
  COFD[         105] =   0.2822883202839590E+01;
  COFD[         106] =  -0.1534915185459436E+00;
  COFD[         107] =   0.6713123822819834E-02;
  COFD[         108] =  -0.1288560752457988E+02;
  COFD[         109] =   0.2706409284361047E+01;
  COFD[         110] =  -0.1381656199059203E+00;
  COFD[         111] =   0.6039972042896966E-02;
  COFD[         112] =  -0.1336434031050532E+02;
  COFD[         113] =   0.2773068323801607E+01;
  COFD[         114] =  -0.1471915976486450E+00;
  COFD[         115] =   0.6445718552318650E-02;
  COFD[         116] =  -0.1417981137873820E+02;
  COFD[         117] =   0.3108160137081816E+01;
  COFD[         118] =  -0.1911804380629321E+00;
  COFD[         119] =   0.8371995333872883E-02;
  COFD[         120] =  -0.1447793881138075E+02;
  COFD[         121] =   0.3898507010768475E+01;
  COFD[         122] =  -0.2933558378767546E+00;
  COFD[         123] =   0.1277297948985323E-01;
  COFD[         124] =  -0.1304821435550623E+02;
  COFD[         125] =   0.2822365020409423E+01;
  COFD[         126] =  -0.1534211775973923E+00;
  COFD[         127] =   0.6709953971249074E-02;
  COFD[         128] =  -0.1306285439260799E+02;
  COFD[         129] =   0.2821999540521922E+01;
  COFD[         130] =  -0.1533715125660936E+00;
  COFD[         131] =   0.6707714606075711E-02;
  COFD[         132] =  -0.1060239225167790E+02;
  COFD[         133] =   0.2155926512014339E+01;
  COFD[         134] =  -0.6501121368580048E-01;
  COFD[         135] =   0.2795359361150588E-02;
  COFD[         136] =  -0.1443057934118857E+02;
  COFD[         137] =   0.3202458275062821E+01;
  COFD[         138] =  -0.2029322289873638E+00;
  COFD[         139] =   0.8859791257831449E-02;
  COFD[         140] =  -0.1444835589212445E+02;
  COFD[         141] =   0.3207617516919998E+01;
  COFD[         142] =  -0.2036261056164307E+00;
  COFD[         143] =   0.8890786328902988E-02;
  COFD[         144] =  -0.1878205175654568E+02;
  COFD[         145] =   0.4857153319943412E+01;
  COFD[         146] =  -0.4000141640745818E+00;
  COFD[         147] =   0.1668590773597161E-01;
  COFD[         148] =  -0.1446569243405954E+02;
  COFD[         149] =   0.3212703400637984E+01;
  COFD[         150] =  -0.2043100955529739E+00;
  COFD[         151] =   0.8921339110627818E-02;
  COFD[         152] =  -0.1421469420975667E+02;
  COFD[         153] =   0.3109400981453699E+01;
  COFD[         154] =  -0.1913493546425900E+00;
  COFD[         155] =   0.8379626477869216E-02;
  COFD[         156] =  -0.1665586079319310E+02;
  COFD[         157] =   0.4002786058866318E+01;
  COFD[         158] =  -0.3040981874501080E+00;
  COFD[         159] =   0.1312476752287860E-01;
  COFD[         160] =  -0.1719172783314806E+02;
  COFD[         161] =   0.4153756740871656E+01;
  COFD[         162] =  -0.3216968524975696E+00;
  COFD[         163] =   0.1380380031784617E-01;
  COFD[         164] =  -0.1299956565222219E+02;
  COFD[         165] =   0.2824331144452877E+01;
  COFD[         166] =  -0.1536882294633172E+00;
  COFD[         167] =   0.6721992209348116E-02;
  COFD[         168] =  -0.1290038354247348E+02;
  COFD[         169] =   0.2706727627372812E+01;
  COFD[         170] =  -0.1382090269907308E+00;
  COFD[         171] =   0.6041930389082575E-02;
  COFD[         172] =  -0.1336450146586354E+02;
  COFD[         173] =   0.2764941911287849E+01;
  COFD[         174] =  -0.1460911037285431E+00;
  COFD[         175] =   0.6396245185837268E-02;
  COFD[         176] =  -0.1417599919891023E+02;
  COFD[         177] =   0.3098479314022822E+01;
  COFD[         178] =  -0.1898625531460105E+00;
  COFD[         179] =   0.8312456268065094E-02;
  COFD[         180] =  -0.1127823399275531E+02;
  COFD[         181] =   0.2707816433593029E+01;
  COFD[         182] =  -0.1386141258544173E+00;
  COFD[         183] =   0.6070036436918267E-02;
  COFD[         184] =  -0.1059938809341299E+02;
  COFD[         185] =   0.2155942578374456E+01;
  COFD[         186] =  -0.6501315584613947E-01;
  COFD[         187] =   0.2795438090048244E-02;
  COFD[         188] =  -0.1060239225167790E+02;
  COFD[         189] =   0.2155926512014339E+01;
  COFD[         190] =  -0.6501121368580048E-01;
  COFD[         191] =   0.2795359361150588E-02;
  COFD[         192] =  -0.9780465322652098E+01;
  COFD[         193] =   0.1955037452782714E+01;
  COFD[         194] =  -0.4069384010078646E-01;
  COFD[         195] =   0.1831517429116509E-02;
  COFD[         196] =  -0.1191485716730014E+02;
  COFD[         197] =   0.2573525435966725E+01;
  COFD[         198] =  -0.1209966797925984E+00;
  COFD[         199] =   0.5292806824954136E-02;
  COFD[         200] =  -0.1191832199726468E+02;
  COFD[         201] =   0.2574628563741546E+01;
  COFD[         202] =  -0.1211478456187282E+00;
  COFD[         203] =   0.5299681334394376E-02;
  COFD[         204] =  -0.1755770485245349E+02;
  COFD[         205] =   0.4765099994814697E+01;
  COFD[         206] =  -0.3967432989876266E+00;
  COFD[         207] =   0.1688051763142629E-01;
  COFD[         208] =  -0.1192160240338415E+02;
  COFD[         209] =   0.2575674668634583E+01;
  COFD[         210] =  -0.1212911957444766E+00;
  COFD[         211] =   0.5306200358730889E-02;
  COFD[         212] =  -0.1187870349337985E+02;
  COFD[         213] =   0.2544526033876783E+01;
  COFD[         214] =  -0.1172898003477687E+00;
  COFD[         215] =   0.5133388723630184E-02;
  COFD[         216] =  -0.1396088659212954E+02;
  COFD[         217] =   0.3362678569206413E+01;
  COFD[         218] =  -0.2260211902233529E+00;
  COFD[         219] =   0.9953713808622193E-02;
  COFD[         220] =  -0.1467854583311947E+02;
  COFD[         221] =   0.3596468340692584E+01;
  COFD[         222] =  -0.2558538870403423E+00;
  COFD[         223] =   0.1122178946707502E-01;
  COFD[         224] =  -0.1056125561830246E+02;
  COFD[         225] =   0.2156060195364817E+01;
  COFD[         226] =  -0.6502865120836085E-01;
  COFD[         227] =   0.2796118937582714E-02;
  COFD[         228] =  -0.1098942473992561E+02;
  COFD[         229] =   0.2261791237155641E+01;
  COFD[         230] =  -0.8129994784319448E-01;
  COFD[         231] =   0.3616450326455660E-02;
  COFD[         232] =  -0.1116880542243876E+02;
  COFD[         233] =   0.2261180518893006E+01;
  COFD[         234] =  -0.8121943300255434E-01;
  COFD[         235] =   0.3612953520751428E-02;
  COFD[         236] =  -0.1186868669178899E+02;
  COFD[         237] =   0.2541192297100982E+01;
  COFD[         238] =  -0.1168329960776270E+00;
  COFD[         239] =   0.5112616487592947E-02;
  COFD[         240] =  -0.1642636080356586E+02;
  COFD[         241] =   0.4526058244629064E+01;
  COFD[         242] =  -0.3728782124089335E+00;
  COFD[         243] =   0.1613117546672936E-01;
  COFD[         244] =  -0.1443514465604553E+02;
  COFD[         245] =   0.3212703400637984E+01;
  COFD[         246] =  -0.2043100955529739E+00;
  COFD[         247] =   0.8921339110627811E-02;
  COFD[         248] =  -0.1443057934118857E+02;
  COFD[         249] =   0.3202458275062821E+01;
  COFD[         250] =  -0.2029322289873638E+00;
  COFD[         251] =   0.8859791257831449E-02;
  COFD[         252] =  -0.1191485716730014E+02;
  COFD[         253] =   0.2573525435966725E+01;
  COFD[         254] =  -0.1209966797925984E+00;
  COFD[         255] =   0.5292806824954136E-02;
  COFD[         256] =  -0.1549200610767626E+02;
  COFD[         257] =   0.3439251591988564E+01;
  COFD[         258] =  -0.2323545954069675E+00;
  COFD[         259] =   0.1007912823695402E-01;
  COFD[         260] =  -0.1550002724603720E+02;
  COFD[         261] =   0.3439385567461067E+01;
  COFD[         262] =  -0.2323724410704331E+00;
  COFD[         263] =   0.1007991807289719E-01;
  COFD[         264] =  -0.2035995576301296E+02;
  COFD[         265] =   0.5195752187472611E+01;
  COFD[         266] =  -0.4301703237047207E+00;
  COFD[         267] =   0.1745457159373258E-01;
  COFD[         268] =  -0.1550831727319405E+02;
  COFD[         269] =   0.3439771059431501E+01;
  COFD[         270] =  -0.2324237888991663E+00;
  COFD[         271] =   0.1008219068775275E-01;
  COFD[         272] =  -0.1525512227936506E+02;
  COFD[         273] =   0.3342973466357380E+01;
  COFD[         274] =  -0.2204616885047087E+00;
  COFD[         275] =   0.9589723244771213E-02;
  COFD[         276] =  -0.1751574793768704E+02;
  COFD[         277] =   0.4133936172799989E+01;
  COFD[         278] =  -0.3179073696597811E+00;
  COFD[         279] =   0.1359390538028741E-01;
  COFD[         280] =  -0.1800123749909511E+02;
  COFD[         281] =   0.4268166153537729E+01;
  COFD[         282] =  -0.3331465709119429E+00;
  COFD[         283] =   0.1416531900068370E-01;
  COFD[         284] =  -0.1442811073014737E+02;
  COFD[         285] =   0.3232109859231266E+01;
  COFD[         286] =  -0.2069193497139531E+00;
  COFD[         287] =   0.9037874467945330E-02;
  COFD[         288] =  -0.1422382512076542E+02;
  COFD[         289] =   0.3088771824088037E+01;
  COFD[         290] =  -0.1886569944625107E+00;
  COFD[         291] =   0.8262345406342235E-02;
  COFD[         292] =  -0.1441524433081045E+02;
  COFD[         293] =   0.3015383266367500E+01;
  COFD[         294] =  -0.1786596762821762E+00;
  COFD[         295] =   0.7810473248205813E-02;
  COFD[         296] =  -0.1524045238644407E+02;
  COFD[         297] =   0.3344366702662799E+01;
  COFD[         298] =  -0.2206501613083914E+00;
  COFD[         299] =   0.9598185922907161E-02;
  COFD[         300] =  -0.1643428330282377E+02;
  COFD[         301] =   0.4529060545715668E+01;
  COFD[         302] =  -0.3732727663838795E+00;
  COFD[         303] =   0.1614842241690972E-01;
  COFD[         304] =  -0.1445309446773090E+02;
  COFD[         305] =   0.3218021790448843E+01;
  COFD[         306] =  -0.2050253330481837E+00;
  COFD[         307] =   0.8953286991448825E-02;
  COFD[         308] =  -0.1444835589212445E+02;
  COFD[         309] =   0.3207617516919998E+01;
  COFD[         310] =  -0.2036261056164307E+00;
  COFD[         311] =   0.8890786328902988E-02;
  COFD[         312] =  -0.1191832199726468E+02;
  COFD[         313] =   0.2574628563741546E+01;
  COFD[         314] =  -0.1211478456187282E+00;
  COFD[         315] =   0.5299681334394376E-02;
  COFD[         316] =  -0.1550002724603720E+02;
  COFD[         317] =   0.3439385567461067E+01;
  COFD[         318] =  -0.2323724410704331E+00;
  COFD[         319] =   0.1007991807289719E-01;
  COFD[         320] =  -0.1550751325209437E+02;
  COFD[         321] =   0.3439251591988575E+01;
  COFD[         322] =  -0.2323545954069691E+00;
  COFD[         323] =   0.1007912823695409E-01;
  COFD[         324] =  -0.1969599282188356E+02;
  COFD[         325] =   0.4976350873560208E+01;
  COFD[         326] =  -0.4064271824554607E+00;
  COFD[         327] =   0.1660711621307244E-01;
  COFD[         328] =  -0.1551528528364382E+02;
  COFD[         329] =   0.3439377629946809E+01;
  COFD[         330] =  -0.2323713837861868E+00;
  COFD[         331] =   0.1007987127828529E-01;
  COFD[         332] =  -0.1526366040509268E+02;
  COFD[         333] =   0.3343427896083488E+01;
  COFD[         334] =  -0.2205232081495746E+00;
  COFD[         335] =   0.9592487262706322E-02;
  COFD[         336] =  -0.1751683577397043E+02;
  COFD[         337] =   0.4130764670922778E+01;
  COFD[         338] =  -0.3174968952024533E+00;
  COFD[         339] =   0.1357624209632402E-01;
  COFD[         340] =  -0.1800286708825678E+02;
  COFD[         341] =   0.4265349872271823E+01;
  COFD[         342] =  -0.3327930702686513E+00;
  COFD[         343] =   0.1415058099824088E-01;
  COFD[         344] =  -0.1444671785222707E+02;
  COFD[         345] =   0.3237789274993721E+01;
  COFD[         346] =  -0.2076830584081396E+00;
  COFD[         347] =   0.9071985084775150E-02;
  COFD[         348] =  -0.1423890668903722E+02;
  COFD[         349] =   0.3093149275869819E+01;
  COFD[         350] =  -0.1892529443234314E+00;
  COFD[         351] =   0.8289269066164192E-02;
  COFD[         352] =  -0.1442342778716563E+02;
  COFD[         353] =   0.3015752514175661E+01;
  COFD[         354] =  -0.1787103636683844E+00;
  COFD[         355] =   0.7812778229298878E-02;
  COFD[         356] =  -0.1525024957896720E+02;
  COFD[         357] =   0.3345452630291720E+01;
  COFD[         358] =  -0.2207969417284804E+00;
  COFD[         359] =   0.9604772071746296E-02;
  COFD[         360] =  -0.1734880949753627E+02;
  COFD[         361] =   0.4593172772581599E+01;
  COFD[         362] =  -0.3376039007204080E+00;
  COFD[         363] =   0.1287220770047757E-01;
  COFD[         364] =  -0.1876832382228059E+02;
  COFD[         365] =   0.4857884811542586E+01;
  COFD[         366] =  -0.4000862896404223E+00;
  COFD[         367] =   0.1668799943730030E-01;
  COFD[         368] =  -0.1878205175654568E+02;
  COFD[         369] =   0.4857153319943412E+01;
  COFD[         370] =  -0.4000141640745818E+00;
  COFD[         371] =   0.1668590773597161E-01;
  COFD[         372] =  -0.1755770485245349E+02;
  COFD[         373] =   0.4765099994814697E+01;
  COFD[         374] =  -0.3967432989876266E+00;
  COFD[         375] =   0.1688051763142629E-01;
  COFD[         376] =  -0.2035995576301296E+02;
  COFD[         377] =   0.5195752187472611E+01;
  COFD[         378] =  -0.4301703237047207E+00;
  COFD[         379] =   0.1745457159373258E-01;
  COFD[         380] =  -0.1969599282188356E+02;
  COFD[         381] =   0.4976350873560208E+01;
  COFD[         382] =  -0.4064271824554607E+00;
  COFD[         383] =   0.1660711621307244E-01;
  COFD[         384] =  -0.1387567274892456E+02;
  COFD[         385] =   0.1812858971426756E+01;
  COFD[         386] =   0.1097216544334093E+00;
  COFD[         387] =  -0.9394611941949941E-02;
  COFD[         388] =  -0.1969031169212750E+02;
  COFD[         389] =   0.4971091362886221E+01;
  COFD[         390] =  -0.4056069022661880E+00;
  COFD[         391] =   0.1656594510101084E-01;
  COFD[         392] =  -0.2025023470652300E+02;
  COFD[         393] =   0.5164679115096308E+01;
  COFD[         394] =  -0.4289835854551929E+00;
  COFD[         395] =   0.1751207656231255E-01;
  COFD[         396] =  -0.2010037767510506E+02;
  COFD[         397] =   0.4898106549869066E+01;
  COFD[         398] =  -0.3743028741194254E+00;
  COFD[         399] =   0.1438132798773566E-01;
  COFD[         400] =  -0.2011880612108384E+02;
  COFD[         401] =   0.4811463838735835E+01;
  COFD[         402] =  -0.3574368957363810E+00;
  COFD[         403] =   0.1344766487934238E-01;
  COFD[         404] =  -0.1871632696644692E+02;
  COFD[         405] =   0.4858563271173600E+01;
  COFD[         406] =  -0.4001100488600431E+00;
  COFD[         407] =   0.1668621701604891E-01;
  COFD[         408] =  -0.1866117368121606E+02;
  COFD[         409] =   0.4778493812881500E+01;
  COFD[         410] =  -0.3921938583237028E+00;
  COFD[         411] =   0.1644103604696690E-01;
  COFD[         412] =  -0.1896062192545326E+02;
  COFD[         413] =   0.4755814945087671E+01;
  COFD[         414] =  -0.3879510962583117E+00;
  COFD[         415] =   0.1620356073022963E-01;
  COFD[         416] =  -0.2025517466685841E+02;
  COFD[         417] =   0.5173702332068737E+01;
  COFD[         418] =  -0.4304299361036623E+00;
  COFD[         419] =   0.1758606643266914E-01;
  COFD[         420] =  -0.1644185676789823E+02;
  COFD[         421] =   0.4531933357015490E+01;
  COFD[         422] =  -0.3736503142227166E+00;
  COFD[         423] =   0.1616492637429652E-01;
  COFD[         424] =  -0.1447055666380567E+02;
  COFD[         425] =   0.3223243595288925E+01;
  COFD[         426] =  -0.2057275599882618E+00;
  COFD[         427] =   0.8984653058524278E-02;
  COFD[         428] =  -0.1446569243405954E+02;
  COFD[         429] =   0.3212703400637984E+01;
  COFD[         430] =  -0.2043100955529739E+00;
  COFD[         431] =   0.8921339110627818E-02;
  COFD[         432] =  -0.1192160240338415E+02;
  COFD[         433] =   0.2575674668634583E+01;
  COFD[         434] =  -0.1212911957444766E+00;
  COFD[         435] =   0.5306200358730889E-02;
  COFD[         436] =  -0.1550831727319405E+02;
  COFD[         437] =   0.3439771059431501E+01;
  COFD[         438] =  -0.2324237888991663E+00;
  COFD[         439] =   0.1008219068775275E-01;
  COFD[         440] =  -0.1551528528364382E+02;
  COFD[         441] =   0.3439377629946809E+01;
  COFD[         442] =  -0.2323713837861868E+00;
  COFD[         443] =   0.1007987127828529E-01;
  COFD[         444] =  -0.1969031169212750E+02;
  COFD[         445] =   0.4971091362886221E+01;
  COFD[         446] =  -0.4056069022661880E+00;
  COFD[         447] =   0.1656594510101084E-01;
  COFD[         448] =  -0.1552255388569028E+02;
  COFD[         449] =   0.3439251591988566E+01;
  COFD[         450] =  -0.2323545954069677E+00;
  COFD[         451] =   0.1007912823695403E-01;
  COFD[         452] =  -0.1527248932147248E+02;
  COFD[         453] =   0.3344142969538372E+01;
  COFD[         454] =  -0.2206199135639200E+00;
  COFD[         455] =   0.9596828436646991E-02;
  COFD[         456] =  -0.1751830948101999E+02;
  COFD[         457] =   0.4127897839132661E+01;
  COFD[         458] =  -0.3171257099248815E+00;
  COFD[         459] =   0.1356026356523252E-01;
  COFD[         460] =  -0.1800478358618446E+02;
  COFD[         461] =   0.4262790016299151E+01;
  COFD[         462] =  -0.3324715275574177E+00;
  COFD[         463] =   0.1413716549763151E-01;
  COFD[         464] =  -0.1446477198029847E+02;
  COFD[         465] =   0.3243341387819648E+01;
  COFD[         466] =  -0.2084296261317005E+00;
  COFD[         467] =   0.9105329427825336E-02;
  COFD[         468] =  -0.1425351882632697E+02;
  COFD[         469] =   0.3097428667095943E+01;
  COFD[         470] =  -0.1898355238768725E+00;
  COFD[         471] =   0.8315588096767848E-02;
  COFD[         472] =  -0.1443174973667854E+02;
  COFD[         473] =   0.3016320413575568E+01;
  COFD[         474] =  -0.1787880785865032E+00;
  COFD[         475] =   0.7816303442184303E-02;
  COFD[         476] =  -0.1526028005455884E+02;
  COFD[         477] =   0.3346771853142201E+01;
  COFD[         478] =  -0.2209752142445045E+00;
  COFD[         479] =   0.9612769758541551E-02;
  COFD[         480] =  -0.1655778410328815E+02;
  COFD[         481] =   0.4581623131684943E+01;
  COFD[         482] =  -0.3829686886254905E+00;
  COFD[         483] =   0.1668786826439446E-01;
  COFD[         484] =  -0.1421921207753527E+02;
  COFD[         485] =   0.3119581660542013E+01;
  COFD[         486] =  -0.1927352195932762E+00;
  COFD[         487] =   0.8442234647456850E-02;
  COFD[         488] =  -0.1421469420975667E+02;
  COFD[         489] =   0.3109400981453699E+01;
  COFD[         490] =  -0.1913493546425900E+00;
  COFD[         491] =   0.8379626477869216E-02;
  COFD[         492] =  -0.1187870349337985E+02;
  COFD[         493] =   0.2544526033876783E+01;
  COFD[         494] =  -0.1172898003477687E+00;
  COFD[         495] =   0.5133388723630184E-02;
  COFD[         496] =  -0.1525512227936506E+02;
  COFD[         497] =   0.3342973466357380E+01;
  COFD[         498] =  -0.2204616885047087E+00;
  COFD[         499] =   0.9589723244771213E-02;
  COFD[         500] =  -0.1526366040509268E+02;
  COFD[         501] =   0.3343427896083488E+01;
  COFD[         502] =  -0.2205232081495746E+00;
  COFD[         503] =   0.9592487262706322E-02;
  COFD[         504] =  -0.2025023470652300E+02;
  COFD[         505] =   0.5164679115096308E+01;
  COFD[         506] =  -0.4289835854551929E+00;
  COFD[         507] =   0.1751207656231255E-01;
  COFD[         508] =  -0.1527248932147248E+02;
  COFD[         509] =   0.3344142969538372E+01;
  COFD[         510] =  -0.2206199135639200E+00;
  COFD[         511] =   0.9596828436646991E-02;
  COFD[         512] =  -0.1497768289771445E+02;
  COFD[         513] =   0.3227480066943741E+01;
  COFD[         514] =  -0.2056129721234479E+00;
  COFD[         515] =   0.8952934338714936E-02;
  COFD[         516] =  -0.1727794685544993E+02;
  COFD[         517] =   0.4041230524213696E+01;
  COFD[         518] =  -0.3067591011045543E+00;
  COFD[         519] =   0.1314646875574740E-01;
  COFD[         520] =  -0.1780793471702511E+02;
  COFD[         521] =   0.4199950217760466E+01;
  COFD[         522] =  -0.3259890442365050E+00;
  COFD[         523] =   0.1392425431224151E-01;
  COFD[         524] =  -0.1421309755343899E+02;
  COFD[         525] =   0.3139101796054804E+01;
  COFD[         526] =  -0.1953917148330065E+00;
  COFD[         527] =   0.8562227807514249E-02;
  COFD[         528] =  -0.1395546834860353E+02;
  COFD[         529] =   0.2971056305378044E+01;
  COFD[         530] =  -0.1731879094770783E+00;
  COFD[         531] =   0.7584632150714469E-02;
  COFD[         532] =  -0.1415896169467945E+02;
  COFD[         533] =   0.2906352555651907E+01;
  COFD[         534] =  -0.1643826505243385E+00;
  COFD[         535] =   0.7187021854047564E-02;
  COFD[         536] =  -0.1496176862944729E+02;
  COFD[         537] =   0.3228126895460024E+01;
  COFD[         538] =  -0.2057002163408370E+00;
  COFD[         539] =   0.8956841574723430E-02;
  COFD[         540] =  -0.1765762436087865E+02;
  COFD[         541] =   0.4866441430895612E+01;
  COFD[         542] =  -0.4030149143618096E+00;
  COFD[         543] =   0.1685813093438678E-01;
  COFD[         544] =  -0.1666875714366874E+02;
  COFD[         545] =   0.4017058659430817E+01;
  COFD[         546] =  -0.3059783073756657E+00;
  COFD[         547] =   0.1320711908953701E-01;
  COFD[         548] =  -0.1665586079319310E+02;
  COFD[         549] =   0.4002786058866318E+01;
  COFD[         550] =  -0.3040981874501080E+00;
  COFD[         551] =   0.1312476752287860E-01;
  COFD[         552] =  -0.1396088659212954E+02;
  COFD[         553] =   0.3362678569206413E+01;
  COFD[         554] =  -0.2260211902233529E+00;
  COFD[         555] =   0.9953713808622193E-02;
  COFD[         556] =  -0.1751574793768704E+02;
  COFD[         557] =   0.4133936172799989E+01;
  COFD[         558] =  -0.3179073696597811E+00;
  COFD[         559] =   0.1359390538028741E-01;
  COFD[         560] =  -0.1751683577397043E+02;
  COFD[         561] =   0.4130764670922778E+01;
  COFD[         562] =  -0.3174968952024533E+00;
  COFD[         563] =   0.1357624209632402E-01;
  COFD[         564] =  -0.2010037767510506E+02;
  COFD[         565] =   0.4898106549869066E+01;
  COFD[         566] =  -0.3743028741194254E+00;
  COFD[         567] =   0.1438132798773566E-01;
  COFD[         568] =  -0.1751830948101999E+02;
  COFD[         569] =   0.4127897839132661E+01;
  COFD[         570] =  -0.3171257099248815E+00;
  COFD[         571] =   0.1356026356523252E-01;
  COFD[         572] =  -0.1727794685544993E+02;
  COFD[         573] =   0.4041230524213696E+01;
  COFD[         574] =  -0.3067591011045543E+00;
  COFD[         575] =   0.1314646875574740E-01;
  COFD[         576] =  -0.1948513279043824E+02;
  COFD[         577] =   0.4762989340707079E+01;
  COFD[         578] =  -0.3907037809895256E+00;
  COFD[         579] =   0.1639711559555987E-01;
  COFD[         580] =  -0.1994609678650828E+02;
  COFD[         581] =   0.4879517700734503E+01;
  COFD[         582] =  -0.4026370932876875E+00;
  COFD[         583] =   0.1678808270297160E-01;
  COFD[         584] =  -0.1668569159498162E+02;
  COFD[         585] =   0.4046494293563119E+01;
  COFD[         586] =  -0.3098570121881317E+00;
  COFD[         587] =   0.1337707028099142E-01;
  COFD[         588] =  -0.1636978212746447E+02;
  COFD[         589] =   0.3866683320271529E+01;
  COFD[         590] =  -0.2875792696211273E+00;
  COFD[         591] =   0.1245724760985546E-01;
  COFD[         592] =  -0.1635720575391811E+02;
  COFD[         593] =   0.3700108154273381E+01;
  COFD[         594] =  -0.2653005108968848E+00;
  COFD[         595] =   0.1146722423828080E-01;
  COFD[         596] =  -0.1727623628623193E+02;
  COFD[         597] =   0.4048916194516065E+01;
  COFD[         598] =  -0.3077539539381744E+00;
  COFD[         599] =   0.1318928678668319E-01;
  COFD[         600] =  -0.1786291016737098E+02;
  COFD[         601] =   0.4848948909820512E+01;
  COFD[         602] =  -0.3945292942419953E+00;
  COFD[         603] =   0.1623802827748926E-01;
  COFD[         604] =  -0.1720382549451406E+02;
  COFD[         605] =   0.4167339155910690E+01;
  COFD[         606] =  -0.3234483507589772E+00;
  COFD[         607] =   0.1387892452356678E-01;
  COFD[         608] =  -0.1719172783314806E+02;
  COFD[         609] =   0.4153756740871656E+01;
  COFD[         610] =  -0.3216968524975696E+00;
  COFD[         611] =   0.1380380031784617E-01;
  COFD[         612] =  -0.1467854583311947E+02;
  COFD[         613] =   0.3596468340692584E+01;
  COFD[         614] =  -0.2558538870403423E+00;
  COFD[         615] =   0.1122178946707502E-01;
  COFD[         616] =  -0.1800123749909511E+02;
  COFD[         617] =   0.4268166153537729E+01;
  COFD[         618] =  -0.3331465709119429E+00;
  COFD[         619] =   0.1416531900068370E-01;
  COFD[         620] =  -0.1800286708825678E+02;
  COFD[         621] =   0.4265349872271823E+01;
  COFD[         622] =  -0.3327930702686513E+00;
  COFD[         623] =   0.1415058099824088E-01;
  COFD[         624] =  -0.2011880612108384E+02;
  COFD[         625] =   0.4811463838735835E+01;
  COFD[         626] =  -0.3574368957363810E+00;
  COFD[         627] =   0.1344766487934238E-01;
  COFD[         628] =  -0.1800478358618446E+02;
  COFD[         629] =   0.4262790016299151E+01;
  COFD[         630] =  -0.3324715275574177E+00;
  COFD[         631] =   0.1413716549763151E-01;
  COFD[         632] =  -0.1780793471702511E+02;
  COFD[         633] =   0.4199950217760466E+01;
  COFD[         634] =  -0.3259890442365050E+00;
  COFD[         635] =   0.1392425431224151E-01;
  COFD[         636] =  -0.1994609678650828E+02;
  COFD[         637] =   0.4879517700734503E+01;
  COFD[         638] =  -0.4026370932876875E+00;
  COFD[         639] =   0.1678808270297160E-01;
  COFD[         640] =  -0.2042982002220235E+02;
  COFD[         641] =   0.5005361629731139E+01;
  COFD[         642] =  -0.4157082406938964E+00;
  COFD[         643] =   0.1722721086482850E-01;
  COFD[         644] =  -0.1722254346565480E+02;
  COFD[         645] =   0.4196189217828696E+01;
  COFD[         646] =  -0.3271714819192909E+00;
  COFD[         647] =   0.1403874343755477E-01;
  COFD[         648] =  -0.1701594967403162E+02;
  COFD[         649] =   0.4069211126205633E+01;
  COFD[         650] =  -0.3126069256643104E+00;
  COFD[         651] =   0.1348840574684995E-01;
  COFD[         652] =  -0.1697793472542246E+02;
  COFD[         653] =   0.3898476651990144E+01;
  COFD[         654] =  -0.2900815498039417E+00;
  COFD[         655] =   0.1250041749182838E-01;
  COFD[         656] =  -0.1780986732405179E+02;
  COFD[         657] =   0.4208912464383648E+01;
  COFD[         658] =  -0.3271471801976215E+00;
  COFD[         659] =   0.1397401213591626E-01;
  COFD[         660] =  -0.1434851661908174E+02;
  COFD[         661] =   0.3864395688187810E+01;
  COFD[         662] =  -0.2888220557172488E+00;
  COFD[         663] =   0.1257261544647088E-01;
  COFD[         664] =  -0.1298157331946940E+02;
  COFD[         665] =   0.2822883202839590E+01;
  COFD[         666] =  -0.1534915185459436E+00;
  COFD[         667] =   0.6713123822819834E-02;
  COFD[         668] =  -0.1299956565222219E+02;
  COFD[         669] =   0.2824331144452877E+01;
  COFD[         670] =  -0.1536882294633172E+00;
  COFD[         671] =   0.6721992209348116E-02;
  COFD[         672] =  -0.1056125561830246E+02;
  COFD[         673] =   0.2156060195364817E+01;
  COFD[         674] =  -0.6502865120836085E-01;
  COFD[         675] =   0.2796118937582714E-02;
  COFD[         676] =  -0.1442811073014737E+02;
  COFD[         677] =   0.3232109859231266E+01;
  COFD[         678] =  -0.2069193497139531E+00;
  COFD[         679] =   0.9037874467945330E-02;
  COFD[         680] =  -0.1444671785222707E+02;
  COFD[         681] =   0.3237789274993721E+01;
  COFD[         682] =  -0.2076830584081396E+00;
  COFD[         683] =   0.9071985084775150E-02;
  COFD[         684] =  -0.1871632696644692E+02;
  COFD[         685] =   0.4858563271173600E+01;
  COFD[         686] =  -0.4001100488600431E+00;
  COFD[         687] =   0.1668621701604891E-01;
  COFD[         688] =  -0.1446477198029847E+02;
  COFD[         689] =   0.3243341387819648E+01;
  COFD[         690] =  -0.2084296261317005E+00;
  COFD[         691] =   0.9105329427825336E-02;
  COFD[         692] =  -0.1421309755343899E+02;
  COFD[         693] =   0.3139101796054804E+01;
  COFD[         694] =  -0.1953917148330065E+00;
  COFD[         695] =   0.8562227807514249E-02;
  COFD[         696] =  -0.1668569159498162E+02;
  COFD[         697] =   0.4046494293563119E+01;
  COFD[         698] =  -0.3098570121881317E+00;
  COFD[         699] =   0.1337707028099142E-01;
  COFD[         700] =  -0.1722254346565480E+02;
  COFD[         701] =   0.4196189217828696E+01;
  COFD[         702] =  -0.3271714819192909E+00;
  COFD[         703] =   0.1403874343755477E-01;
  COFD[         704] =  -0.1292646230231139E+02;
  COFD[         705] =   0.2821999540521940E+01;
  COFD[         706] =  -0.1533715125660963E+00;
  COFD[         707] =   0.6707714606075836E-02;
  COFD[         708] =  -0.1284085524690963E+02;
  COFD[         709] =   0.2707949300899787E+01;
  COFD[         710] =  -0.1383738279427983E+00;
  COFD[         711] =   0.6049324358308415E-02;
  COFD[         712] =  -0.1335047052208153E+02;
  COFD[         713] =   0.2788888190516315E+01;
  COFD[         714] =  -0.1493333595043647E+00;
  COFD[         715] =   0.6541990186326518E-02;
  COFD[         716] =  -0.1417184929578018E+02;
  COFD[         717] =   0.3126696914084966E+01;
  COFD[         718] =  -0.1937032186843793E+00;
  COFD[         719] =   0.8485951819761462E-02;
  COFD[         720] =  -0.1473052091803629E+02;
  COFD[         721] =   0.3938367747699902E+01;
  COFD[         722] =  -0.3007580570513714E+00;
  COFD[         723] =   0.1318674711212272E-01;
  COFD[         724] =  -0.1288560752457988E+02;
  COFD[         725] =   0.2706409284361047E+01;
  COFD[         726] =  -0.1381656199059203E+00;
  COFD[         727] =   0.6039972042896966E-02;
  COFD[         728] =  -0.1290038354247348E+02;
  COFD[         729] =   0.2706727627372812E+01;
  COFD[         730] =  -0.1382090269907308E+00;
  COFD[         731] =   0.6041930389082575E-02;
  COFD[         732] =  -0.1098942473992561E+02;
  COFD[         733] =   0.2261791237155641E+01;
  COFD[         734] =  -0.8129994784319448E-01;
  COFD[         735] =   0.3616450326455660E-02;
  COFD[         736] =  -0.1422382512076542E+02;
  COFD[         737] =   0.3088771824088037E+01;
  COFD[         738] =  -0.1886569944625107E+00;
  COFD[         739] =   0.8262345406342235E-02;
  COFD[         740] =  -0.1423890668903722E+02;
  COFD[         741] =   0.3093149275869819E+01;
  COFD[         742] =  -0.1892529443234314E+00;
  COFD[         743] =   0.8289269066164192E-02;
  COFD[         744] =  -0.1866117368121606E+02;
  COFD[         745] =   0.4778493812881500E+01;
  COFD[         746] =  -0.3921938583237028E+00;
  COFD[         747] =   0.1644103604696690E-01;
  COFD[         748] =  -0.1425351882632697E+02;
  COFD[         749] =   0.3097428667095943E+01;
  COFD[         750] =  -0.1898355238768725E+00;
  COFD[         751] =   0.8315588096767848E-02;
  COFD[         752] =  -0.1395546834860353E+02;
  COFD[         753] =   0.2971056305378044E+01;
  COFD[         754] =  -0.1731879094770783E+00;
  COFD[         755] =   0.7584632150714469E-02;
  COFD[         756] =  -0.1636978212746447E+02;
  COFD[         757] =   0.3866683320271529E+01;
  COFD[         758] =  -0.2875792696211273E+00;
  COFD[         759] =   0.1245724760985546E-01;
  COFD[         760] =  -0.1701594967403162E+02;
  COFD[         761] =   0.4069211126205633E+01;
  COFD[         762] =  -0.3126069256643104E+00;
  COFD[         763] =   0.1348840574684995E-01;
  COFD[         764] =  -0.1284085524690963E+02;
  COFD[         765] =   0.2707949300899787E+01;
  COFD[         766] =  -0.1383738279427983E+00;
  COFD[         767] =   0.6049324358308415E-02;
  COFD[         768] =  -0.1267750532012987E+02;
  COFD[         769] =   0.2569994454759337E+01;
  COFD[         770] =  -0.1199417264901710E+00;
  COFD[         771] =   0.5227731201117593E-02;
  COFD[         772] =  -0.1306696965803475E+02;
  COFD[         773] =   0.2610536352821046E+01;
  COFD[         774] =  -0.1255043916260873E+00;
  COFD[         775] =   0.5480787778210062E-02;
  COFD[         776] =  -0.1392335517962079E+02;
  COFD[         777] =   0.2962245288209091E+01;
  COFD[         778] =  -0.1719893475330243E+00;
  COFD[         779] =   0.7530527438410521E-02;
  COFD[         780] =  -0.1555514842196157E+02;
  COFD[         781] =   0.4198163984729581E+01;
  COFD[         782] =  -0.3357679563943392E+00;
  COFD[         783] =   0.1475425817654628E-01;
  COFD[         784] =  -0.1336434031050532E+02;
  COFD[         785] =   0.2773068323801607E+01;
  COFD[         786] =  -0.1471915976486450E+00;
  COFD[         787] =   0.6445718552318650E-02;
  COFD[         788] =  -0.1336450146586354E+02;
  COFD[         789] =   0.2764941911287849E+01;
  COFD[         790] =  -0.1460911037285431E+00;
  COFD[         791] =   0.6396245185837268E-02;
  COFD[         792] =  -0.1116880542243876E+02;
  COFD[         793] =   0.2261180518893006E+01;
  COFD[         794] =  -0.8121943300255434E-01;
  COFD[         795] =   0.3612953520751428E-02;
  COFD[         796] =  -0.1441524433081045E+02;
  COFD[         797] =   0.3015383266367500E+01;
  COFD[         798] =  -0.1786596762821762E+00;
  COFD[         799] =   0.7810473248205813E-02;
  COFD[         800] =  -0.1442342778716563E+02;
  COFD[         801] =   0.3015752514175661E+01;
  COFD[         802] =  -0.1787103636683844E+00;
  COFD[         803] =   0.7812778229298878E-02;
  COFD[         804] =  -0.1896062192545326E+02;
  COFD[         805] =   0.4755814945087671E+01;
  COFD[         806] =  -0.3879510962583117E+00;
  COFD[         807] =   0.1620356073022963E-01;
  COFD[         808] =  -0.1443174973667854E+02;
  COFD[         809] =   0.3016320413575568E+01;
  COFD[         810] =  -0.1787880785865032E+00;
  COFD[         811] =   0.7816303442184303E-02;
  COFD[         812] =  -0.1415896169467945E+02;
  COFD[         813] =   0.2906352555651907E+01;
  COFD[         814] =  -0.1643826505243385E+00;
  COFD[         815] =   0.7187021854047564E-02;
  COFD[         816] =  -0.1635720575391811E+02;
  COFD[         817] =   0.3700108154273381E+01;
  COFD[         818] =  -0.2653005108968848E+00;
  COFD[         819] =   0.1146722423828080E-01;
  COFD[         820] =  -0.1697793472542246E+02;
  COFD[         821] =   0.3898476651990144E+01;
  COFD[         822] =  -0.2900815498039417E+00;
  COFD[         823] =   0.1250041749182838E-01;
  COFD[         824] =  -0.1335047052208153E+02;
  COFD[         825] =   0.2788888190516315E+01;
  COFD[         826] =  -0.1493333595043647E+00;
  COFD[         827] =   0.6541990186326518E-02;
  COFD[         828] =  -0.1306696965803475E+02;
  COFD[         829] =   0.2610536352821046E+01;
  COFD[         830] =  -0.1255043916260873E+00;
  COFD[         831] =   0.5480787778210062E-02;
  COFD[         832] =  -0.1332407043105473E+02;
  COFD[         833] =   0.2569994454759347E+01;
  COFD[         834] =  -0.1199417264901723E+00;
  COFD[         835] =   0.5227731201117647E-02;
  COFD[         836] =  -0.1414287301380644E+02;
  COFD[         837] =   0.2906808913725623E+01;
  COFD[         838] =  -0.1644441595237770E+00;
  COFD[         839] =   0.7189777153290596E-02;
  COFD[         840] =  -0.1653528873346560E+02;
  COFD[         841] =   0.4572908957761951E+01;
  COFD[         842] =  -0.3818057701602784E+00;
  COFD[         843] =   0.1663628005999293E-01;
  COFD[         844] =  -0.1417981137873820E+02;
  COFD[         845] =   0.3108160137081816E+01;
  COFD[         846] =  -0.1911804380629321E+00;
  COFD[         847] =   0.8371995333872883E-02;
  COFD[         848] =  -0.1417599919891023E+02;
  COFD[         849] =   0.3098479314022822E+01;
  COFD[         850] =  -0.1898625531460105E+00;
  COFD[         851] =   0.8312456268065094E-02;
  COFD[         852] =  -0.1186868669178899E+02;
  COFD[         853] =   0.2541192297100982E+01;
  COFD[         854] =  -0.1168329960776270E+00;
  COFD[         855] =   0.5112616487592947E-02;
  COFD[         856] =  -0.1524045238644407E+02;
  COFD[         857] =   0.3344366702662799E+01;
  COFD[         858] =  -0.2206501613083914E+00;
  COFD[         859] =   0.9598185922907161E-02;
  COFD[         860] =  -0.1525024957896720E+02;
  COFD[         861] =   0.3345452630291720E+01;
  COFD[         862] =  -0.2207969417284804E+00;
  COFD[         863] =   0.9604772071746296E-02;
  COFD[         864] =  -0.2025517466685841E+02;
  COFD[         865] =   0.5173702332068737E+01;
  COFD[         866] =  -0.4304299361036623E+00;
  COFD[         867] =   0.1758606643266914E-01;
  COFD[         868] =  -0.1526028005455884E+02;
  COFD[         869] =   0.3346771853142201E+01;
  COFD[         870] =  -0.2209752142445045E+00;
  COFD[         871] =   0.9612769758541551E-02;
  COFD[         872] =  -0.1496176862944729E+02;
  COFD[         873] =   0.3228126895460024E+01;
  COFD[         874] =  -0.2057002163408370E+00;
  COFD[         875] =   0.8956841574723430E-02;
  COFD[         876] =  -0.1727623628623193E+02;
  COFD[         877] =   0.4048916194516065E+01;
  COFD[         878] =  -0.3077539539381744E+00;
  COFD[         879] =   0.1318928678668319E-01;
  COFD[         880] =  -0.1780986732405179E+02;
  COFD[         881] =   0.4208912464383648E+01;
  COFD[         882] =  -0.3271471801976215E+00;
  COFD[         883] =   0.1397401213591626E-01;
  COFD[         884] =  -0.1417184929578018E+02;
  COFD[         885] =   0.3126696914084966E+01;
  COFD[         886] =  -0.1937032186843793E+00;
  COFD[         887] =   0.8485951819761462E-02;
  COFD[         888] =  -0.1392335517962079E+02;
  COFD[         889] =   0.2962245288209091E+01;
  COFD[         890] =  -0.1719893475330243E+00;
  COFD[         891] =   0.7530527438410521E-02;
  COFD[         892] =  -0.1414287301380644E+02;
  COFD[         893] =   0.2906808913725623E+01;
  COFD[         894] =  -0.1644441595237770E+00;
  COFD[         895] =   0.7189777153290596E-02;
  COFD[         896] =  -0.1494332402346468E+02;
  COFD[         897] =   0.3227480066943751E+01;
  COFD[         898] =  -0.2056129721234494E+00;
  COFD[         899] =   0.8952934338715007E-02;
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
  COFTD[           4] =   0.2976868391871377E+00;
  COFTD[           5] =   0.2854588501200347E-03;
  COFTD[           6] =  -0.1285543991944700E-06;
  COFTD[           7] =   0.1757513172747589E-10;
  COFTD[           8] =   0.2999265878772748E+00;
  COFTD[           9] =   0.2876065973546457E-03;
  COFTD[          10] =  -0.1295216221593593E-06;
  COFTD[          11] =   0.1770736423857070E-10;
  COFTD[          12] =   0.1585931146622087E+00;
  COFTD[          13] =   0.3775348885882738E-04;
  COFTD[          14] =  -0.1768016190115858E-07;
  COFTD[          15] =   0.2636753851315786E-11;
  COFTD[          16] =   0.2536890768339781E+00;
  COFTD[          17] =   0.3891614006021053E-03;
  COFTD[          18] =  -0.1706164447255217E-06;
  COFTD[          19] =   0.2277134100195441E-10;
  COFTD[          20] =   0.2541780967931188E+00;
  COFTD[          21] =   0.3899115617623597E-03;
  COFTD[          22] =  -0.1709453309663866E-06;
  COFTD[          23] =   0.2281523583725962E-10;
  COFTD[          24] =  -0.1112921453625776E+00;
  COFTD[          25] =   0.6841397958433722E-03;
  COFTD[          26] =  -0.2506086153883979E-06;
  COFTD[          27] =   0.2969326279265050E-10;
  COFTD[          28] =   0.2546389683061027E+00;
  COFTD[          29] =   0.3906185429447143E-03;
  COFTD[          30] =  -0.1712552862076685E-06;
  COFTD[          31] =   0.2285660404479628E-10;
  COFTD[          32] =   0.2738502817645918E+00;
  COFTD[          33] =   0.3598280636877769E-03;
  COFTD[          34] =  -0.1591730965616415E-06;
  COFTD[          35] =   0.2140486083562457E-10;
  COFTD[          36] =   0.1108921212448803E+00;
  COFTD[          37] =   0.5679286832661998E-03;
  COFTD[          38] =  -0.2344616821568033E-06;
  COFTD[          39] =   0.2994319284275358E-10;
  COFTD[          40] =   0.7464590665177903E-01;
  COFTD[          41] =   0.6023631134912639E-03;
  COFTD[          42] =  -0.2449185295932349E-06;
  COFTD[          43] =   0.3096529767229981E-10;
  COFTD[          44] =   0.2952265988608936E+00;
  COFTD[          45] =   0.2830996683152013E-03;
  COFTD[          46] =  -0.1274919581477715E-06;
  COFTD[          47] =   0.1742988161184074E-10;
  COFTD[          48] =   0.3135317345848753E+00;
  COFTD[          49] =   0.2503523445647176E-03;
  COFTD[          50] =  -0.1138565566730279E-06;
  COFTD[          51] =   0.1572601886907115E-10;
  COFTD[          52] =   0.3378441547112856E+00;
  COFTD[          53] =   0.2697655991392449E-03;
  COFTD[          54] =  -0.1226854187454628E-06;
  COFTD[          55] =   0.1694547302788842E-10;
  COFTD[          56] =   0.2725430730178022E+00;
  COFTD[          57] =   0.3581104448883287E-03;
  COFTD[          58] =  -0.1584132928369996E-06;
  COFTD[          59] =   0.2130268594966920E-10;
  COFTD[          60] =  -0.1585931146622087E+00;
  COFTD[          61] =  -0.3775348885882738E-04;
  COFTD[          62] =   0.1768016190115858E-07;
  COFTD[          63] =  -0.2636753851315786E-11;
  COFTD[          64] =   0.4204727183103451E+00;
  COFTD[          65] =  -0.3287526079537624E-05;
  COFTD[          66] =   0.6103160956596783E-08;
  COFTD[          67] =  -0.8093596164801657E-12;
  COFTD[          68] =   0.4268965221228311E+00;
  COFTD[          69] =  -0.3337751508308041E-05;
  COFTD[          70] =   0.6196402460537293E-08;
  COFTD[          71] =  -0.8217246693447408E-12;
  COFTD[          72] =   0.0000000000000000E+00;
  COFTD[          73] =   0.0000000000000000E+00;
  COFTD[          74] =   0.0000000000000000E+00;
  COFTD[          75] =   0.0000000000000000E+00;
  COFTD[          76] =   0.4534800422581371E+00;
  COFTD[          77] =   0.4085412451125554E-04;
  COFTD[          78] =  -0.1740693371737612E-07;
  COFTD[          79] =   0.2810982215515255E-11;
  COFTD[          80] =   0.4552350965244072E+00;
  COFTD[          81] =   0.4101223776616561E-04;
  COFTD[          82] =  -0.1747430187128919E-07;
  COFTD[          83] =   0.2821861252892914E-11;
  COFTD[          84] =   0.9225454658923754E-01;
  COFTD[          85] =   0.4742311371310592E-03;
  COFTD[          86] =  -0.1957520914900890E-06;
  COFTD[          87] =   0.2499719240888261E-10;
  COFTD[          88] =   0.4568919543649144E+00;
  COFTD[          89] =   0.4116150448179982E-04;
  COFTD[          90] =  -0.1753790073324839E-07;
  COFTD[          91] =   0.2832131601065486E-11;
  COFTD[          92] =   0.4585173035952478E+00;
  COFTD[          93] =   0.2448303357627405E-04;
  COFTD[          94] =  -0.8803935139397828E-08;
  COFTD[          95] =   0.1529314166116477E-11;
  COFTD[          96] =   0.3873786628706396E+00;
  COFTD[          97] =   0.1817636775569728E-03;
  COFTD[          98] =  -0.8458974202644200E-07;
  COFTD[          99] =   0.1206323436896060E-10;
  COFTD[         100] =   0.3599424148527454E+00;
  COFTD[         101] =   0.2205256508391225E-03;
  COFTD[         102] =  -0.1015802405503637E-06;
  COFTD[         103] =   0.1425033065176923E-10;
  COFTD[         104] =   0.4134627185516987E+00;
  COFTD[         105] =  -0.3232717393930960E-05;
  COFTD[         106] =   0.6001410819264168E-08;
  COFTD[         107] =  -0.7958662066366476E-12;
  COFTD[         108] =   0.4111024136789401E+00;
  COFTD[         109] =  -0.1436862257324850E-04;
  COFTD[         110] =   0.1257740607483724E-07;
  COFTD[         111] =  -0.1888957757243028E-11;
  COFTD[         112] =   0.4779741554107228E+00;
  COFTD[         113] =  -0.1670588644178477E-04;
  COFTD[         114] =   0.1462330272420397E-07;
  COFTD[         115] =  -0.2196224003028831E-11;
  COFTD[         116] =   0.4541344777639836E+00;
  COFTD[         117] =   0.2424900779110863E-04;
  COFTD[         118] =  -0.8719781032141305E-08;
  COFTD[         119] =   0.1514695922532610E-11;
};




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!
! ************************************************************************
! ************************************************************************
! ******************* Master mechanism C/H/O/N       ********************
! ************************************************************************
! ************************************************************************
!
! SKR/GLA04:  Ø. Skreiberg, P. Kilpinen and P. Glarborg; “Ammonia Chemistry 
!             under Fuel-Rich Conditions in a Flow Reactor”, Combust. Flame
!             136, 501-508 (2004)
!
! RAS/GLA08 : C.L. Rasmussen, J. Hansen, P. Marshall, and P. Glarborg; 
!             “Experimental Measurements and Kinetic Modeling of CO/H2/O2/NOx 
!             Conversion at High Pressure”, Int. J. Chem. Kin. 40, 454-480 (2008)
!
 ELEMENTS
 O  H  N 
 END  
!
 SPECIES    
H  O  OH  H2  O2  HO2  H2O  H2O2            
!
NO NO2 N2O  
!
NH N NNH N2
 END
!
 THERMO ALL
   300.00   1000.00   5000.00  
H                 L 6/94H   1    0    0    0G   200.00   6000.00  1000.        1
 0.25000000E+01 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 0.25473660E+05-0.44668285E+00 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.25473660E+05-0.44668285E+00 0.26219035E+05    4
HO2               L 5/89H   1O   2    0    0G   200.00   6000.00  1000.        1
 4.17226590E+00 1.88120980E-03-3.46292970E-07 1.94685160E-11 1.76091530E-16    2
 3.02010736E+01 2.95697380E+00 4.30178800E+00-4.74902010E-03 2.11579530E-05    3
-2.42759610E-08 9.29206700E-12 2.63190983E+02 3.71587740E+00                   4
H2  REF ELEMENT   RUS 78H   2    0    0    0G   200.00   6000.00  1000.        1
 0.29328305E+01 0.82659802E-03-0.14640057E-06 0.15409851E-10-0.68879615E-15    2
-0.81305582E+03-0.10243164E+01 0.23443029E+01 0.79804248E-02-0.19477917E-04    3
 0.20156967E-07-0.73760289E-11-0.91792413E+03 0.68300218E+00 0.00000000E+00    4
H2O               L 5/89H   2O   1    0    0G   200.00   6000.00  1000.        1
 0.26770389E+01 0.29731816E-02-0.77376889E-06 0.94433514E-10-0.42689991E-14    2
-0.29885894E+05 0.68825500E+01 0.41986352E+01-0.20364017E-02 0.65203416E-05    3
-0.54879269E-08 0.17719680E-11-0.30293726E+05-0.84900901E+00-0.29084817E+05    4
H2O2              T 8/03H   2O   2    0    0G   200.00   6000.00  1000.        1
 4.57977305E+00 4.05326003E-03-1.29844730E-06 1.98211400E-10-1.13968792E-14    2
-1.79847939E+04 6.64969660E-01 4.31515149E+00-8.47390622E-04 1.76404323E-05    3
-2.26762944E-08 9.08950158E-12-1.76843601E+04 3.27373216E+00                   4
N                 120186N   1               G  0300.00   5000.00  1000.00      1
 0.02450268E+02 0.01066146E-02-0.07465337E-06 0.01879652E-09-0.01025984E-13    2
 0.05611604E+06 0.04448758E+02 0.02503071E+02-0.02180018E-03 0.05420529E-06    3
-0.05647560E-09 0.02099904E-12 0.05609890E+06 0.04167566E+02                   4
NH                      N   1H   1    0    0G   200.00   6000.00  1000.0       1
 0.27836929E+01 0.13298429E-02-0.42478047E-06 0.78348504E-10-0.55044470E-14    2
 0.42134514E+05 0.57407798E+01 0.34929084E+01 0.31179197E-03-0.14890484E-05    3
 0.24816442E-08-0.10356967E-11 0.41894294E+05 0.18483277E+01 0.42940822E+05    4
NNH               120186N   2H   1          G  0250.000  4000.00  1000.00      1
 0.04415342E+02 0.01614388E-01-0.01632894E-05-0.08559846E-09 0.01614791E-12    2
 0.02788029E+06 0.09042888E+01 0.03501344E+02 0.02053587E-01 0.07170410E-05    3
 0.04921348E-08-0.09671170E-11 0.02833347E+06 0.06391837E+02                   4
NO                RUS 89N   1O   1    0    0G   200.00   6000.00  1000.        1
 3.26071234E+00 1.19101135E-03-4.29122646E-07 6.94481463E-11-4.03295681E-15    2
 9.89456954E+03 6.36900469E+00 4.21859896E+00-4.63988124E-03 1.10443049E-05    3
-9.34055507E-09 2.80554874E-12 9.81823786E+03 2.28060952E+00                   4
NO2               L 7/88N   1O   2    0    0G   200.00   6000.00  1000.        1
 4.88475400E+00 2.17239550E-03-8.28069090E-07 1.57475100E-10-1.05108950E-14    2
 2.29397777E+03-1.17416951E-01 3.94403120E+00-1.58542900E-03 1.66578120E-05    3
-2.04754260E-08 7.83505640E-12 2.87409757E+03 6.31199190E+00                   4
N2        BUR0302 G 8/02N  2.   0.   0.   0.G   200.00   6000.00  1000.        1
 2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2
-9.23948688E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3
 2.43530612E-09-1.40881235E-12-1.04697628E+03 2.96747038E+00 0.00000000E+00    4
N2O               121286N   2O   1          G  0300.00   5000.00  1000.00      1
 0.04718977E+02 0.02873714E-01-0.01197496E-04 0.02250552E-08-0.01575337E-12    2
 0.08165811E+05-0.01657250E+02 0.02543058E+02 0.09492193E-01-0.09792775E-04    3
 0.06263845E-07-0.01901826E-10 0.08765100E+05 0.09511222E+02                   4
O                 L 1/90O   1    0    0    0G   200.00   6000.00  1000.        1
 2.54363697E+00-2.73162486E-05-4.19029520E-09 4.95481845E-12-4.79553694E-16    2
 2.92260120E+04 4.92229457E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00 2.99687009E+04    4
O2 REF ELEMENT    RUS 89O   2    0    0    0G   200.00   6000.00  1000.        1
 3.66096083E+00 6.56365523E-04-1.41149485E-07 2.05797658E-11-1.29913248E-15    2
-1.21597725E+03 3.41536184E+00 3.78245636E+00-2.99673415E-03 9.84730200E-06    3
-9.68129508E-09 3.24372836E-12-1.06394356E+03 3.65767573E+00 0.00000000E+00    4
OH HYDROXYL RADI  IU3/03O   1H   1    0    0G   200.00   6000.00  1000.        1
 2.83853033E+00 1.10741289E-03-2.94000209E-07 4.20698729E-11-2.42289890E-15    2
 3.70056220E+03 5.84513094E+00 3.99198424E+00-2.40106655E-03 4.61664033E-06    3
-3.87916306E-09 1.36319502E-12 3.37165248E+03-1.03814059E-01                   4
!
 END
!
 REACTIONS
!
! *****************************************************************************
!    H2/O2 subset                                                             *
! *****************************************************************************
!
   H+O2=O+OH                            3.6E15  -0.410   16600 ! RAS/GLA08  HES98
   H+H+M=H2+M                           7.0E17  -1.000       0 ! RAS/GLA08  COH/WES83
     N2/0/ H2O/0/ H2/0/                                        !
   H+H+N2=H2+N2                         5.4E18  -1.300       0 ! RAS/GLA08  COH/WES83
   H+H+H2=H2+H2                         1.0E17  -0.600       0 ! RAS/GLA08  COH/WES83
   H+H+H2O=H2+H2O                       1.0E19  -1.000       0 ! RAS/GLA08  COH/WES83
   H+O+M=OH+M                           6.2E16  -0.600       0 ! RAS/GLA08  MIL/BOW89
     H2O/5/                                                    ! 
   H+O2(+M)=HO2(+M)                     1.5E12   0.600       0 ! RAS/GLA08  MUE/DRY98
     LOW  / 3.5E16 -0.41 -1116 /                               !
     TROE / 0.5 1.0E-30 1.0E30 /                               !
     N2/0/ H2O/11/ H2/2/ O2/0.78/                        !
   H+O2(+N2)=HO2(+N2)                   1.5E12   0.600       0 ! RAS/GLA08  LI/DRY04
     LOW  / 6.37E20 -1.720 520 /                               !
     TROE / 0.8 1.0E-30 1.0E30 /                               !
   O+O+M=O2+M                           1.9E13   0.000   -1788 ! RAS/GLA08  NBS86
     N2/1.5/ O2/1.5/ H2O/10/                                   !
   O+H2=OH+H                            3.8E12   0.000    7948 ! RAS/GLA08  CEC05
     DUPLICATE                                                 !
   O+H2=OH+H                            8.8E14   0.000   19175 ! RAS/GLA08  CEC05
     DUPLICATE                                                 !
   OH+OH=O+H2O                          4.3E03   2.700   -1822 ! RAS/GLA08  SRI/MIC06
   OH+H+M=H2O+M                         4.5E22  -2.000       0 ! RAS/GLA08  CON/WES04
     H2/0.73/ H2O/12/ !HE/0.38/                       !   
   OH+H2=H+H2O                          2.1E08   1.520    3449 ! RAS/GLA08  MIC92
   H2+O2=HO2+H                          7.4E05   2.433   53502 ! RAS/GLA08  MIC/WAG00 
   HO2+H=OH+OH                          8.4E13   0.000     400 ! RAS/GLA08  RAS/GLA08 
   HO2+H=H2O+O                          1.4E12   0.000       0 ! RAS/GLA08  CEC05
   HO2+O=OH+O2                          1.6E13   0.000    -445 ! RAS/GLA08  CEC05
   HO2+OH=H2O+O2                        3.6E21  -2.100    9000 ! RAS/GLA08  RAS/GLA08 
     DUPLICATE                                                 !
   HO2+OH=H2O+O2                        2.0E15  -0.600       0 ! 
     DUPLICATE                                                 !
   HO2+OH=H2O+O2                       -2.2E96 -24.000   49000 !
     DUPLICATE                                                 !
   HO2+HO2=H2O2+O2                      1.9E11   0.000   -1408 ! RAS/GLA08  KAP/TROE02
     DUPLICATE                                                 !
   HO2+HO2=H2O2+O2                      1.0E14   0.000   11034 !
     DUPLICATE                                                 !
   H2O2(+M)=OH+OH(+M)                   4.0E11   0.000   37137 ! RAS/GLA08  KAP/TRO02
     LOW  /2.291E16 0.0 43638/                                 !
     TROE /0.5 1E-30 1E30 1E30/                                !  (Fc=0.5)
     H2O/12/ H2/2.5/                                   !
   H2O2+H=H2O+OH                        1.0E13   0.000    3580 ! RAS/GLA08  CEC05
   H2O2+H=HO2+H2                        1.7E12   0.000    3760 ! RAS/GLA08  CEC05
   H2O2+O=HO2+OH                        9.6E06   2.000    3970 ! RAS/GLA08  NBS86
   H2O2+OH=H2O+HO2                      1.9E12   0.000     427 ! RAS/GLA08  HIP/TRO95
     DUPLICATE                                                 !
   H2O2+OH=H2O+HO2                      1.6E18   0.000   29410 ! RAS/GLA08  HIP/TRO95
     DUPLICATE                                                 !
!
! ************************************************************************
!    N subset (oxid.)                                                    *
! ************************************************************************
!
   NO+O(+M)=NO2(+M)                     1.3E15  -0.750       0 ! RAS/GLA08 ALL/DRY97,NBS91
     LOW  /4.72E24 -2.87 1550/                                 ! RAS/GLA08 ALL/DRY97 (Fc=0.95-1E-04*T)
     TROE /0.880 1E03 1E04 1E30/                               ! RAS/GLA08  (  1bar) ¤
   NO+HO2=NO2+OH                        2.1E12   0.000    -497 ! RAS/GLA08 CEC05
!
   NO2+H=NO+OH                          1.3E14   0.000     362 ! RAS/GLA08 KO/FON91
   NO2+O=NO+O2                          1.1E14  -0.520       0 ! RAS/GLA08 BEM/CLY74
   NO2+NO2=NO+NO+O2                     4.5E12   0.000   27599 ! RAS/GLA08 PAR/LIN98
!
   N2O(+M)=N2+O(+M)                     1.3E12   0.000   62570 ! SKR/GLA04 JOH/GLA92,ROH/HAN96
     LOW/4.0E14 0 56600/
     N2/1.7/ O2/1.4/ H2O/12/   !
   N2O+H=N2+OH                          3.3E10   0.000    4729 ! SKR/GLA04 MAR/FON87
     DUP
   N2O+H=N2+OH                          4.4E14   0.000   19254 !
     DUP
   N2O+O=NO+NO                          9.2E13   0.000   27679 ! SKR/GLA04 MEA/AND00
   N2O+O=N2+O2                          3.7E12   0.000   15936 ! SKR/GLA04 MEA/AND00
!
   NH+H=N+H2                            3.0E13   0.000       0 ! SKR/GLA04 DAV/HAN90,rv
   NH+O=NO+H                            9.2E13   0.000       0 ! SKR/GLA04 CEC94
   NH+OH=N+H2O                          5.0E11   0.500    2000 ! SKR/GLA04 JAM est
   NH+O2=NO+OH                          1.3E06   1.500     100 ! SKR/GLA04 MIL/MEL92
   NH+NO=N2O+H                          2.9E14  -0.400       0 ! SKR/GLA04 MIL/MEL92
     DUP
   NH+NO=N2O+H                         -2.2E13  -0.230       0 !
     DUP
   NH+NO=N2+OH                          2.2E13  -0.230       0 ! SKR/GLA04 MIL/MEL92
   NH+NO2=N2O+OH                        1.0E13   0.000       0 ! SKR/GLA04 HAR/PHI86
   N+OH=NO+H                            3.8E13   0.000       0 ! SKR/GLA04 FLO/HAN77,HOW/SMI80
   N+O2=NO+O                            6.4E09   1.000    6280 ! SKR/GLA04 BAU/DRY73
   N+NO=N2+O                            2.1E13   0.000       0 !           CEC05
!   
   NNH=N2+H                             6.5E07   0.000       0 ! SKR/GLA04 MIL/GLA99
   NNH+H=N2+H2                          1.0E14   0.000       0 ! SKR/GLA04 JAM est
   NNH+O=N2O+H                          1.0E14   0.000       0 ! SKR/GLA04 JAM est
   NNH+O=N2+OH                          8.0E13   0.000       0 ! SKR/GLA04 JAM est
   NNH+O=NH+NO                          5.0E13   0.000       0 ! SKR/GLA04 MIL/MEL92
   NNH+OH=N2+H2O                        5.0E13   0.000       0 ! SKR/GLA04 JAM est
   NNH+O2=N2+HO2                        2.0E14   0.000       0 ! SKR/GLA04 MIL/GLA99
   NNH+O2=N2+H+O2                       5.0E13   0.000       0 ! SKR/GLA04 MIL/GLA99
!
 END


\\
\\
\\  This is the therm file
\\
\\
 THERMO 
   300.00   1000.00   5000.00  
H                 L 6/94H   1    0    0    0G   200.00   6000.00  1000.        1
 0.25000000E+01 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 0.25473660E+05-0.44668285E+00 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.25473660E+05-0.44668285E+00 0.26219035E+05    4
HO2               L 5/89H   1O   2    0    0G   200.00   6000.00  1000.        1
 4.17226590E+00 1.88120980E-03-3.46292970E-07 1.94685160E-11 1.76091530E-16    2
 3.02010736E+01 2.95697380E+00 4.30178800E+00-4.74902010E-03 2.11579530E-05    3
-2.42759610E-08 9.29206700E-12 2.63190983E+02 3.71587740E+00                   4
H2  REF ELEMENT   RUS 78H   2    0    0    0G   200.00   6000.00  1000.        1
 0.29328305E+01 0.82659802E-03-0.14640057E-06 0.15409851E-10-0.68879615E-15    2
-0.81305582E+03-0.10243164E+01 0.23443029E+01 0.79804248E-02-0.19477917E-04    3
 0.20156967E-07-0.73760289E-11-0.91792413E+03 0.68300218E+00 0.00000000E+00    4
H2O               L 5/89H   2O   1    0    0G   200.00   6000.00  1000.        1
 0.26770389E+01 0.29731816E-02-0.77376889E-06 0.94433514E-10-0.42689991E-14    2
-0.29885894E+05 0.68825500E+01 0.41986352E+01-0.20364017E-02 0.65203416E-05    3
-0.54879269E-08 0.17719680E-11-0.30293726E+05-0.84900901E+00-0.29084817E+05    4
H2O2              T 8/03H   2O   2    0    0G   200.00   6000.00  1000.        1
 4.57977305E+00 4.05326003E-03-1.29844730E-06 1.98211400E-10-1.13968792E-14    2
-1.79847939E+04 6.64969660E-01 4.31515149E+00-8.47390622E-04 1.76404323E-05    3
-2.26762944E-08 9.08950158E-12-1.76843601E+04 3.27373216E+00                   4
N                 120186N   1               G  0300.00   5000.00  1000.00      1
 0.02450268E+02 0.01066146E-02-0.07465337E-06 0.01879652E-09-0.01025984E-13    2
 0.05611604E+06 0.04448758E+02 0.02503071E+02-0.02180018E-03 0.05420529E-06    3
-0.05647560E-09 0.02099904E-12 0.05609890E+06 0.04167566E+02                   4
NH                      N   1H   1    0    0G   200.00   6000.00  1000.0       1
 0.27836929E+01 0.13298429E-02-0.42478047E-06 0.78348504E-10-0.55044470E-14    2
 0.42134514E+05 0.57407798E+01 0.34929084E+01 0.31179197E-03-0.14890484E-05    3
 0.24816442E-08-0.10356967E-11 0.41894294E+05 0.18483277E+01 0.42940822E+05    4
NNH               120186N   2H   1          G  0250.000  4000.00  1000.00      1
 0.04415342E+02 0.01614388E-01-0.01632894E-05-0.08559846E-09 0.01614791E-12    2
 0.02788029E+06 0.09042888E+01 0.03501344E+02 0.02053587E-01 0.07170410E-05    3
 0.04921348E-08-0.09671170E-11 0.02833347E+06 0.06391837E+02                   4
NO                RUS 89N   1O   1    0    0G   200.00   6000.00  1000.        1
 3.26071234E+00 1.19101135E-03-4.29122646E-07 6.94481463E-11-4.03295681E-15    2
 9.89456954E+03 6.36900469E+00 4.21859896E+00-4.63988124E-03 1.10443049E-05    3
-9.34055507E-09 2.80554874E-12 9.81823786E+03 2.28060952E+00                   4
NO2               L 7/88N   1O   2    0    0G   200.00   6000.00  1000.        1
 4.88475400E+00 2.17239550E-03-8.28069090E-07 1.57475100E-10-1.05108950E-14    2
 2.29397777E+03-1.17416951E-01 3.94403120E+00-1.58542900E-03 1.66578120E-05    3
-2.04754260E-08 7.83505640E-12 2.87409757E+03 6.31199190E+00                   4
N2        BUR0302 G 8/02N  2.   0.   0.   0.G   200.00   6000.00  1000.        1
 2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2
-9.23948688E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3
 2.43530612E-09-1.40881235E-12-1.04697628E+03 2.96747038E+00 0.00000000E+00    4
N2O               121286N   2O   1          G  0300.00   5000.00  1000.00      1
 0.04718977E+02 0.02873714E-01-0.01197496E-04 0.02250552E-08-0.01575337E-12    2
 0.08165811E+05-0.01657250E+02 0.02543058E+02 0.09492193E-01-0.09792775E-04    3
 0.06263845E-07-0.01901826E-10 0.08765100E+05 0.09511222E+02                   4
O                 L 1/90O   1    0    0    0G   200.00   6000.00  1000.        1
 2.54363697E+00-2.73162486E-05-4.19029520E-09 4.95481845E-12-4.79553694E-16    2
 2.92260120E+04 4.92229457E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00 2.99687009E+04    4
O2 REF ELEMENT    RUS 89O   2    0    0    0G   200.00   6000.00  1000.        1
 3.66096083E+00 6.56365523E-04-1.41149485E-07 2.05797658E-11-1.29913248E-15    2
-1.21597725E+03 3.41536184E+00 3.78245636E+00-2.99673415E-03 9.84730200E-06    3
-9.68129508E-09 3.24372836E-12-1.06394356E+03 3.65767573E+00 0.00000000E+00    4
OH HYDROXYL RADI  IU3/03O   1H   1    0    0G   200.00   6000.00  1000.        1
 2.83853033E+00 1.10741289E-03-2.94000209E-07 4.20698729E-11-2.42289890E-15    2
 3.70056220E+03 5.84513094E+00 3.99198424E+00-2.40106655E-03 4.61664033E-06    3
-3.87916306E-09 1.36319502E-12 3.37165248E+03-1.03814059E-01                   4
!
 END

\\
\\
\\  This is the tran file
\\
\\
!
!  Transport Properties obtained from GRI-Mech 3.0
!
!  Additions (at top) for "Glarborg" mechanism.
!  Some of these represent changes to the GRI data.
!  Entered by Joe Grcar as recommended by Chris Pope.
!
C2H5CHO            2   424.600     4.820     0.000     0.000     1.000 ! NMM
C2H5CO             2   424.600     4.820     2.700     0.000     1.000 ! NMM, CJP
CH2CN              2   422.220     5.329     3.500     0.000     1.000 ! CJP
CH2HCO             2   436.000     3.970     0.000     0.000     2.000 ! CH2CHO
CH3CN              2   422.220     5.329     3.500     0.000     1.000 ! CJP
CH3HCO             2   436.000     3.970     0.000     0.000     2.000 ! CH3CHO
HONO               2   350.000     3.950     1.639     0.000     1.000 ! CJP
NO3                2   400.000     4.200     0.200     0.000     1.000 ! CJP
OCHCHO             2   440.200     4.010     0.000     0.000     2.000 ! NMM CHOCHO
NCN                1   422.220     5.329     2.000     0.000     1.000 ! CJP
NO                 1   139.320     3.339     0.200     1.760     4.000 ! CJP
HNO                2   170.000     3.430     1.620     0.000     1.000 ! CJP
NO2                2   333.590     3.852     0.400     0.000     1.000 ! CJP
!
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
