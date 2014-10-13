
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
void vcomp_wdot(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
                double * restrict k_f_s, double * restrict Kc_s,
                double * restrict tc, double * restrict invT, double * restrict T);

/* Inverse molecular weights */
static const double imw[3] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 44.009950,  /*CO2 */
    1.0 / 28.013400};  /*N2 */



static double fwd_A[0], fwd_beta[0], fwd_Ea[0];
static double low_A[0], low_beta[0], low_Ea[0];
static double rev_A[0], rev_beta[0], rev_Ea[0];
static double troe_a[0],troe_Ts[0], troe_Tss[0], troe_Tsss[0];
static double sri_a[0], sri_b[0], sri_c[0], sri_d[0], sri_e[0];
static double activation_units[0], prefactor_units[0], phase_units[0];
static int is_PD[0], troe_len[0], sri_len[0], nTB[0], *TBid[0];
static double *TB[0];

static double fwd_A_DEF[0], fwd_beta_DEF[0], fwd_Ea_DEF[0];
static double low_A_DEF[0], low_beta_DEF[0], low_Ea_DEF[0];
static double rev_A_DEF[0], rev_beta_DEF[0], rev_Ea_DEF[0];
static double troe_a_DEF[0],troe_Ts_DEF[0], troe_Tss_DEF[0], troe_Tsss_DEF[0];
static double sri_a_DEF[0], sri_b_DEF[0], sri_c_DEF[0], sri_d_DEF[0], sri_e_DEF[0];
static double activation_units_DEF[0], prefactor_units_DEF[0], phase_units_DEF[0];
static int is_PD_DEF[0], troe_len_DEF[0], sri_len_DEF[0], nTB_DEF[0], *TBid_DEF[0];
static double *TB_DEF[0];
static int rxn_map[0] = {};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<0; ++i) {
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
  if (reaction_id<0 || reaction_id>=0) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=3) {
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
    for (int i=0; i<0; i++) {
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
    for (int i=0; i<0; i++) {
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
  for (int i=0; i<0; ++i) {
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
    SetAllDefaults();
}



/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 3;
    *ii = 0;
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
    for (i=0; i<lenkname*3; i++) {
        kname[i] = ' ';
    }

    /* H2  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 1*lenkname + 0 ] = 'C';
    kname[ 1*lenkname + 1 ] = 'O';
    kname[ 1*lenkname + 2 ] = '2';
    kname[ 1*lenkname + 3 ] = ' ';

    /* N2  */
    kname[ 2*lenkname + 0 ] = 'N';
    kname[ 2*lenkname + 1 ] = '2';
    kname[ 2*lenkname + 2 ] = ' ';

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
    XW += x[1]*44.009950; /*CO2 */
    XW += x[2]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*CO2 */
    YOW += y[2]*imw[2]; /*N2 */
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

    for (int n=0; n<3; n++) {
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
    W += c[1]*44.009950; /*CO2 */
    W += c[2]*28.013400; /*N2 */

    for (id = 0; id < 3; ++id) {
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
    XW += x[1]*44.009950; /*CO2 */
    XW += x[2]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double YOW = 0;
    double tmp[3];

    for (int i = 0; i < 3; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 3; i++)
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
    W += c[1]*44.009950; /*CO2 */
    W += c[2]*28.013400; /*N2 */

    for (id = 0; id < 3; ++id) {
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
    double tmp[3];

    for (int i = 0; i < 3; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 3; i++)
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
    XW += x[1]*44.009950; /*CO2 */
    XW += x[2]*28.013400; /*N2 */
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
    W += c[1]*44.009950; /*CO2 */
    W += c[2]*28.013400; /*N2 */

    for (id = 0; id < 3; ++id) {
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
    double tmp[3];

    for (int i = 0; i < 3; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 3; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 3; i++)
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

    for (int n=0; n<3; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<3; n++) {
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
    for (int i = 0; i < 3; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 3; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 3; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int i = 0; i < 3; i++)
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
    XW += x[1]*44.009950; /*CO2 */
    XW += x[2]*28.013400; /*N2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*44.009950*XWinv; 
    y[2] = x[2]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 3; ++id) {
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
    XW += x[1]*44.009950; /*CO2 */
    XW += x[2]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 3; ++id) {
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
    CW += c[1]*44.009950; /*CO2 */
    CW += c[2]*28.013400; /*N2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*44.009950*CWinv; 
    y[2] = c[2]*28.013400*CWinv; 

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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    cvms[1] *= 1.889234139098090e+06; /*CO2 */
    cvms[2] *= 2.968047434442088e+06; /*N2 */
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
    cpms[1] *= 1.889234139098090e+06; /*CO2 */
    cpms[2] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 3; i++)
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
    for (int i = 0; i < 3; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tc[5], h[3];

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
    }

    for (int n=0; n<3; n++) {
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
    for (int i = 0; i < 3; i++)
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
    for (int i = 0; i < 3; i++)
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
    sms[1] *= 1.889234139098090e+06; /*CO2 */
    sms[2] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[3]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 3; ++id) {
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
    double cpor[3], tresult[3]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 3; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 3; i++)
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
    double cvor[3]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 3; ++id) {
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
    double cvor[3]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*CO2 */
    result += cvor[2]*y[2]*imw[2]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[3]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 3; ++id) {
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
    double hml[3], tmp[3]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 3; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 3; ++id) {
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
    double uml[3]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 3; ++id) {
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
    double ums[3]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*CO2 */
    result += y[2]*ums[2]*imw[2]; /*N2 */

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
    double sor[3]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 3; ++id) {
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
    double sor[3]; /* temporary storage */
    double x[3]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*CO2 */
    YOW += y[2]*imw[2]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(44.009950*YOW); 
    x[2] = y[2]/(28.013400*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
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
    double gort[3]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 3; ++id) {
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
    double gort[3]; /* temporary storage */
    double x[3]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*CO2 */
    YOW += y[2]*imw[2]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(44.009950*YOW); 
    x[2] = y[2]/(28.013400*YOW); 
    gibbs(gort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(gort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(gort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(gort[2]+log((x[2]+1e-100))+logPratio);
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
    double aort[3]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 3; ++id) {
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
    double aort[3]; /* temporary storage */
    double x[3]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*CO2 */
    YOW += y[2]*imw[2]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(44.009950*YOW); 
    x[2] = y[2]/(28.013400*YOW); 
    helmholtz(aort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(aort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(aort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(aort[2]+log((x[2]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 3; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*CO2 */
    YOW += y[2]*imw[2]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 3; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[3*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<3; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<3*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*44.009950; /*CO2 */
    XW += x[2]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 3; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 3; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 3; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*CO2 */
    YOW += y[2]*imw[2]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 3; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*44.009950; /*CO2 */
    XW += x[2]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 3; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
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
    for (id = 0; id < 3 * kd; ++ id) {
         nuki[id] = 0; 
    }
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 3; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 1 ] = 2; /*H */

    /*CO2 */
    ncf[ 1 * kd + 2 ] = 1; /*C */
    ncf[ 1 * kd + 0 ] = 2; /*O */

    /*N2 */
    ncf[ 2 * kd + 3 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)
{
    for (int i=0; i<0; ++i) {
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
    double gort[3]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[3]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[3]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[3]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[3]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[0];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[0];
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

    double qdot, q_f[0], q_r[0];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 3; ++i) {
        wdot[i] = 0.0;
    }

    return;
}

void comp_k_f(double * restrict tc, double invT, double * restrict k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<0; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double * restrict tc, double invT, double * restrict Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[3];
    gibbs(g_RT, tc);


#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<0; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;


    return;
}

void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)
{

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 3; ++i) {
        mixture += sc[i];
    }

    double Corr[0];
    for (int i = 0; i < 0; ++i) {
        Corr[i] = 1.0;
    }

    for (int i=0; i<0; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[0*npt], Kc_s[0*npt], mixture[npt], g_RT[3*npt];
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

    for (int n=0; n<3; n++) {
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
    }
}

void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[3];
        tg[0] = tc[0*npt+i];
        tg[1] = tc[1*npt+i];
        tg[2] = tc[2*npt+i];
        tg[3] = tc[3*npt+i];
        tg[4] = tc[4*npt+i];

        gibbs(g, tg);

        g_RT[0*npt+i] = g[0];
        g_RT[1*npt+i] = g[1];
        g_RT[2*npt+i] = g[2];
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
    }
}

/*compute the reaction Jacobian */
void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)
{
    double c[3];

    for (int k=0; k<3; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    for (int k=0; k<3; k++) {
        J[12+k] *= 1.e-6;
    }

    /* dTdot/d[X] */
    for (int k=0; k<3; k++) {
        J[k*4+3] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<16; i++) {
        J[i] = 0.0;
    }

    double wdot[3];
    for (int k=0; k<3; k++) {
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
    for (int k = 0; k < 3; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[3];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[3];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[3];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    double c_R[3], dcRdT[3], e_RT[3];
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
    for (int k = 0; k < 3; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[12+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 3; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 3; ++m) {
            dehmixdc += eh_RT[m]*J[k*4+m];
        }
        J[k*4+3] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[15] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
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
        /*species 1: CO2 */
        species[1] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
        /*species 2: N2 */
        species[2] =
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

    double q_f[0], q_r[0];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 0; ++i) {
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
        /*species 1: CO2 */
        species[1] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 2: N2 */
        species[2] =
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
        /*species 1: CO2 */
        species[1] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 2: N2 */
        species[2] =
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
    wt[1] = 44.009950; /*CO2 */
    wt[2] = 28.013400; /*N2 */

    return;
}


/*save atomic weights into array */
void atomicWeight(double * restrict awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 12.011150; /*C */
    awt[3] = 14.006700; /*N */

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
  *LENIMC =           13;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC) {
  *LENRMC =          264;}
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
  *KK =            3;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
void egtransetNLITE(int* NLITE) {
  *NLITE =            1;}
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
  WT[           1] =   0.4400995063781738E+02;
  WT[           2] =   0.2801339912414551E+02;
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
  EPS[           1] =   0.2440000000000000E+03;
  EPS[           2] =   0.9753000000000000E+02;
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
  SIG[           1] =   0.3763000000000000E+01;
  SIG[           2] =   0.3621000000000000E+01;
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
  POL[           1] =   0.2650000000000000E+01;
  POL[           2] =   0.1760000000000000E+01;
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
  ZROT[           1] =   0.2100000000000000E+01;
  ZROT[           2] =   0.4000000000000000E+01;
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
  NLIN[           2] =            1;
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
  COFLAM[           4] =  -0.1239079949009015E+02;
  COFLAM[           5] =   0.6341842149327256E+01;
  COFLAM[           6] =  -0.6370214552008658E+00;
  COFLAM[           7] =   0.2372277065304921E-01;
  COFLAM[           8] =   0.1154286034252058E+02;
  COFLAM[           9] =  -0.2911509179774523E+01;
  COFLAM[          10] =   0.5546559973428081E+00;
  COFLAM[          11] =  -0.2750092774174382E-01;
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
  COFETA[           4] =  -0.2364265995035323E+02;
  COFETA[           5] =   0.4983971927098391E+01;
  COFETA[           6] =  -0.5507551343971877E+00;
  COFETA[           7] =   0.2334595082438979E-01;
  COFETA[           8] =  -0.1626172843728175E+02;
  COFETA[           9] =   0.2251740453876138E+01;
  COFETA[          10] =  -0.2138340893699383E+00;
  COFETA[          11] =   0.9477823154448534E-02;
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
  COFD[           4] =  -0.1530635401605573E+02;
  COFD[           5] =   0.3873896336632167E+01;
  COFD[           6] =  -0.2956113676790514E+00;
  COFD[           7] =   0.1311329899072374E-01;
  COFD[           8] =  -0.1217511489930491E+02;
  COFD[           9] =   0.2679286720576265E+01;
  COFD[          10] =  -0.1373997446113497E+00;
  COFD[          11] =   0.6125379451516252E-02;
  COFD[          12] =  -0.1530635401605573E+02;
  COFD[          13] =   0.3873896336632167E+01;
  COFD[          14] =  -0.2956113676790514E+00;
  COFD[          15] =   0.1311329899072374E-01;
  COFD[          16] =  -0.2069000320825284E+02;
  COFD[          17] =   0.5101158792846945E+01;
  COFD[          18] =  -0.4264320198306029E+00;
  COFD[          19] =   0.1763100216755543E-01;
  COFD[          20] =  -0.1811180157901078E+02;
  COFD[          21] =   0.4336077093491970E+01;
  COFD[          22] =  -0.3444079746717316E+00;
  COFD[          23] =   0.1476188131006474E-01;
  COFD[          24] =  -0.1217511489930491E+02;
  COFD[          25] =   0.2679286720576265E+01;
  COFD[          26] =  -0.1373997446113497E+00;
  COFD[          27] =   0.6125379451516252E-02;
  COFD[          28] =  -0.1811180157901078E+02;
  COFD[          29] =   0.4336077093491970E+01;
  COFD[          30] =  -0.3444079746717316E+00;
  COFD[          31] =   0.1476188131006474E-01;
  COFD[          32] =  -0.1521415383275665E+02;
  COFD[          33] =   0.3348053449783496E+01;
  COFD[          34] =  -0.2233657260417595E+00;
  COFD[          35] =   0.9817787279109837E-02;
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
  COFTD[           4] =   0.3257425690224397E+00;
  COFTD[           5] =   0.3036334352955890E-03;
  COFTD[           6] =  -0.1552903446031957E-06;
  COFTD[           7] =   0.2414664576246333E-10;
  COFTD[           8] =   0.4452620885099525E+00;
  COFTD[           9] =   0.4946972054970773E-04;
  COFTD[          10] =  -0.2630235132337235E-07;
  COFTD[          11] =   0.4903063330400376E-11;
};




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!     
!  simple 3 species gas mixture  to fluctest     
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!     
ELEMENTS                                                                        
O  H  C  N                                                                
END                                                                             
SPECIES                                                                         
H2      CO2     N2
END                                                                             
REACTIONS                                                                       
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
