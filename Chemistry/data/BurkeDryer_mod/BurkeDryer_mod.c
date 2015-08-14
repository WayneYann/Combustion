
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
static const double imw[11] = {
    1.0 / 1.007970,  /*H */
    1.0 / 2.015940,  /*H2 */
    1.0 / 15.999400,  /*O */
    1.0 / 17.007370,  /*OH */
    1.0 / 18.015340,  /*H2O */
    1.0 / 31.998800,  /*O2 */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 28.013400,  /*N2 */
    1.0 / 39.948000,  /*AR */
    1.0 / 4.002600};  /*HE */



static double fwd_A[29], fwd_beta[29], fwd_Ea[29];
static double low_A[29], low_beta[29], low_Ea[29];
static double rev_A[29], rev_beta[29], rev_Ea[29];
static double troe_a[29],troe_Ts[29], troe_Tss[29], troe_Tsss[29];
static double sri_a[29], sri_b[29], sri_c[29], sri_d[29], sri_e[29];
static double activation_units[29], prefactor_units[29], phase_units[29];
static int is_PD[29], troe_len[29], sri_len[29], nTB[29], *TBid[29];
static double *TB[29];

static double fwd_A_DEF[29], fwd_beta_DEF[29], fwd_Ea_DEF[29];
static double low_A_DEF[29], low_beta_DEF[29], low_Ea_DEF[29];
static double rev_A_DEF[29], rev_beta_DEF[29], rev_Ea_DEF[29];
static double troe_a_DEF[29],troe_Ts_DEF[29], troe_Tss_DEF[29], troe_Tsss_DEF[29];
static double sri_a_DEF[29], sri_b_DEF[29], sri_c_DEF[29], sri_d_DEF[29], sri_e_DEF[29];
static double activation_units_DEF[29], prefactor_units_DEF[29], phase_units_DEF[29];
static int is_PD_DEF[29], troe_len_DEF[29], sri_len_DEF[29], nTB_DEF[29], *TBid_DEF[29];
static double *TB_DEF[29];
static int rxn_map[29] = {7,8,9,10,11,2,12,13,3,14,15,4,5,16,0,17,18,19,20,21,22,1,23,24,25,26,27,28,6};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<29; ++i) {
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
  if (reaction_id<0 || reaction_id>=29) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=11) {
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
    for (int i=0; i<29; i++) {
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
    for (int i=0; i<29; i++) {
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
  for (int i=0; i<29; ++i) {
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
    fwd_A[7]     = 104000000000000;
    fwd_beta[7]  = 0;
    fwd_Ea[7]    = 15286;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = 1e-12;
    is_PD[7] = 0;
    nTB[7] = 0;

    // (1):  O + H2 <=> H + OH
    fwd_A[8]     = 3818000000000;
    fwd_beta[8]  = 0;
    fwd_Ea[8]    = 7948;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = 1e-12;
    is_PD[8] = 0;
    nTB[8] = 0;

    // (2):  O + H2 <=> H + OH
    fwd_A[9]     = 879200000000000;
    fwd_beta[9]  = 0;
    fwd_Ea[9]    = 19170;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = 1e-12;
    is_PD[9] = 0;
    nTB[9] = 0;

    // (3):  H2 + OH <=> H2O + H
    fwd_A[10]     = 216000000;
    fwd_beta[10]  = 1.51;
    fwd_Ea[10]    = 3430;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = 1e-12;
    is_PD[10] = 0;
    nTB[10] = 0;

    // (4):  OH + OH <=> O + H2O
    fwd_A[11]     = 33400;
    fwd_beta[11]  = 2.4199999999999999;
    fwd_Ea[11]    = -1930;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = 1e-12;
    is_PD[11] = 0;
    nTB[11] = 0;

    // (5):  H2 + M <=> H + H + M
    fwd_A[2]     = 4.577e+19;
    fwd_beta[2]  = -1.3999999999999999;
    fwd_Ea[2]    = 104380;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = 1e-6;
    is_PD[2] = 0;
    nTB[2] = 4;
    TB[2] = (double *) malloc(4 * sizeof(double));
    TBid[2] = (int *) malloc(4 * sizeof(int));
    TBid[2][0] = 1; TB[2][0] = 2.5; // H2
    TBid[2][1] = 4; TB[2][1] = 12; // H2O
    TBid[2][2] = 9; TB[2][2] = 0; // AR
    TBid[2][3] = 10; TB[2][3] = 0; // HE

    // (6):  H2 + AR <=> H + H + AR
    fwd_A[12]     = 5.84e+18;
    fwd_beta[12]  = -1.1000000000000001;
    fwd_Ea[12]    = 104380;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = 1e-12;
    is_PD[12] = 0;
    nTB[12] = 0;

    // (7):  H2 + HE <=> H + H + HE
    fwd_A[13]     = 5.84e+18;
    fwd_beta[13]  = -1.1000000000000001;
    fwd_Ea[13]    = 104380;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = 1e-12;
    is_PD[13] = 0;
    nTB[13] = 0;

    // (8):  O + O + M <=> O2 + M
    fwd_A[3]     = 6165000000000000;
    fwd_beta[3]  = -0.5;
    fwd_Ea[3]    = 0;
    prefactor_units[3]  = 1.0000000000000002e-12;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = 1e-12;
    is_PD[3] = 0;
    nTB[3] = 4;
    TB[3] = (double *) malloc(4 * sizeof(double));
    TBid[3] = (int *) malloc(4 * sizeof(int));
    TBid[3][0] = 1; TB[3][0] = 2.5; // H2
    TBid[3][1] = 4; TB[3][1] = 12; // H2O
    TBid[3][2] = 9; TB[3][2] = 0; // AR
    TBid[3][3] = 10; TB[3][3] = 0; // HE

    // (9):  O + O + AR <=> O2 + AR
    fwd_A[14]     = 18860000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = -1788;
    prefactor_units[14]  = 1.0000000000000002e-12;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = 1e-18;
    is_PD[14] = 0;
    nTB[14] = 0;

    // (10):  O + O + HE <=> O2 + HE
    fwd_A[15]     = 18860000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = -1788;
    prefactor_units[15]  = 1.0000000000000002e-12;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = 1e-18;
    is_PD[15] = 0;
    nTB[15] = 0;

    // (11):  O + H + M <=> OH + M
    fwd_A[4]     = 4.714e+18;
    fwd_beta[4]  = -1;
    fwd_Ea[4]    = 0;
    prefactor_units[4]  = 1.0000000000000002e-12;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = 1e-12;
    is_PD[4] = 0;
    nTB[4] = 4;
    TB[4] = (double *) malloc(4 * sizeof(double));
    TBid[4] = (int *) malloc(4 * sizeof(int));
    TBid[4][0] = 1; TB[4][0] = 2.5; // H2
    TBid[4][1] = 4; TB[4][1] = 12; // H2O
    TBid[4][2] = 9; TB[4][2] = 0.75; // AR
    TBid[4][3] = 10; TB[4][3] = 0.75; // HE

    // (12):  H2O + M <=> H + OH + M
    fwd_A[5]     = 6.0640000000000002e+27;
    fwd_beta[5]  = -3.3220000000000001;
    fwd_Ea[5]    = 120790;
    prefactor_units[5]  = 1.0000000000000002e-06;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = 1e-6;
    is_PD[5] = 0;
    nTB[5] = 5;
    TB[5] = (double *) malloc(5 * sizeof(double));
    TBid[5] = (int *) malloc(5 * sizeof(int));
    TBid[5][0] = 1; TB[5][0] = 3; // H2
    TBid[5][1] = 4; TB[5][1] = 0; // H2O
    TBid[5][2] = 10; TB[5][2] = 1.1000000000000001; // HE
    TBid[5][3] = 8; TB[5][3] = 2; // N2
    TBid[5][4] = 5; TB[5][4] = 1.5; // O2

    // (13):  H2O + H2O <=> H + OH + H2O
    fwd_A[16]     = 1.006e+26;
    fwd_beta[16]  = -2.4399999999999999;
    fwd_Ea[16]    = 120180;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = 1e-12;
    is_PD[16] = 0;
    nTB[16] = 0;

    // (14):  H + O2 (+M) <=> HO2 (+M)
    fwd_A[0]     = 4650840000000;
    fwd_beta[0]  = 0.44;
    fwd_Ea[0]    = 0;
    low_A[0]     = 6.366e+20;
    low_beta[0]  = -1.72;
    low_Ea[0]    = 524.79999999999995;
    troe_a[0]    = 0.5;
    troe_Tsss[0] = 1.0000000000000001e-30;
    troe_Ts[0]   = 1e+30;
    troe_len[0]  = 3;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = 1e-12;
    is_PD[0] = 1;
    nTB[0] = 5;
    TB[0] = (double *) malloc(5 * sizeof(double));
    TBid[0] = (int *) malloc(5 * sizeof(int));
    TBid[0][0] = 1; TB[0][0] = 2; // H2
    TBid[0][1] = 4; TB[0][1] = 14; // H2O
    TBid[0][2] = 5; TB[0][2] = 0.78000000000000003; // O2
    TBid[0][3] = 9; TB[0][3] = 0.67000000000000004; // AR
    TBid[0][4] = 10; TB[0][4] = 0.80000000000000004; // HE

    // (15):  HO2 + H <=> H2 + O2
    fwd_A[17]     = 2750000;
    fwd_beta[17]  = 2.0899999999999999;
    fwd_Ea[17]    = -1451;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = 1e-12;
    is_PD[17] = 0;
    nTB[17] = 0;

    // (16):  HO2 + H <=> OH + OH
    fwd_A[18]     = 70790000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 295;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = 1e-12;
    is_PD[18] = 0;
    nTB[18] = 0;

    // (17):  HO2 + O <=> O2 + OH
    fwd_A[19]     = 28500000000;
    fwd_beta[19]  = 1;
    fwd_Ea[19]    = -723.92999999999995;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = 1e-12;
    is_PD[19] = 0;
    nTB[19] = 0;

    // (18):  HO2 + OH <=> H2O + O2
    fwd_A[20]     = 28900000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = -497;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = 1e-12;
    is_PD[20] = 0;
    nTB[20] = 0;

    // (19):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[21]     = 420000000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = 11982;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = 1e-12;
    is_PD[21] = 0;
    nTB[21] = 0;

    // (20):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[22]     = 130000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = -1629.3;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = 1e-12;
    is_PD[22] = 0;
    nTB[22] = 0;

    // (21):  H2O2 (+M) <=> OH + OH (+M)
    fwd_A[1]     = 2000000000000;
    fwd_beta[1]  = 0.90000000000000002;
    fwd_Ea[1]    = 48749;
    low_A[1]     = 2.49e+24;
    low_beta[1]  = -2.2999999999999998;
    low_Ea[1]    = 48749;
    troe_a[1]    = 0.42999999999999999;
    troe_Tsss[1] = 1.0000000000000001e-30;
    troe_Ts[1]   = 1e+30;
    troe_len[1]  = 3;
    prefactor_units[1]  = 1;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = 1e-6;
    is_PD[1] = 1;
    nTB[1] = 6;
    TB[1] = (double *) malloc(6 * sizeof(double));
    TBid[1] = (int *) malloc(6 * sizeof(int));
    TBid[1][0] = 4; TB[1][0] = 7.5; // H2O
    TBid[1][1] = 8; TB[1][1] = 1.5; // N2
    TBid[1][2] = 5; TB[1][2] = 1.2; // O2
    TBid[1][3] = 10; TB[1][3] = 0.65000000000000002; // HE
    TBid[1][4] = 7; TB[1][4] = 7.7000000000000002; // H2O2
    TBid[1][5] = 1; TB[1][5] = 3.7000000000000002; // H2

    // (22):  H2O2 + H <=> H2O + OH
    fwd_A[23]     = 24100000000000;
    fwd_beta[23]  = 0;
    fwd_Ea[23]    = 3970;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = 1e-12;
    is_PD[23] = 0;
    nTB[23] = 0;

    // (23):  H2O2 + H <=> HO2 + H2
    fwd_A[24]     = 48200000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 7950;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = 1e-12;
    is_PD[24] = 0;
    nTB[24] = 0;

    // (24):  H2O2 + O <=> OH + HO2
    fwd_A[25]     = 9550000;
    fwd_beta[25]  = 2;
    fwd_Ea[25]    = 3970;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = 1e-12;
    is_PD[25] = 0;
    nTB[25] = 0;

    // (25):  H2O2 + OH <=> HO2 + H2O
    fwd_A[26]     = 1740000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 318;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = 1e-12;
    is_PD[26] = 0;
    nTB[26] = 0;

    // (26):  H2O2 + OH <=> HO2 + H2O
    fwd_A[27]     = 75900000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = 7270;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = 1e-12;
    is_PD[27] = 0;
    nTB[27] = 0;

    // (27):  HO2 + H <=> O + H2O
    fwd_A[28]     = 3970000000000;
    fwd_beta[28]  = 0;
    fwd_Ea[28]    = 671;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = 1e-12;
    is_PD[28] = 0;
    nTB[28] = 0;

    // (28):  O + OH + M <=> HO2 + M
    fwd_A[6]     = 8000000000000000;
    fwd_beta[6]  = 0;
    fwd_Ea[6]    = 0;
    prefactor_units[6]  = 1.0000000000000002e-12;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = 1e-12;
    is_PD[6] = 0;
    nTB[6] = 4;
    TB[6] = (double *) malloc(4 * sizeof(double));
    TBid[6] = (int *) malloc(4 * sizeof(int));
    TBid[6][0] = 1; TB[6][0] = 2; // H2
    TBid[6][1] = 4; TB[6][1] = 12; // H2O
    TBid[6][2] = 9; TB[6][2] = 0.69999999999999996; // AR
    TBid[6][3] = 10; TB[6][3] = 0.69999999999999996; // HE

    SetAllDefaults();
}



/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 6;
    *kk = 11;
    *ii = 29;
    *nfit = -1; /*Why do you need this anyway ?  */
}

char *strtok_r (char *s, const char *delim, char **save_ptr)
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
    char *p; /*String Tokens */
    char cstr[1000];
    char *saveptr;
    /* Strip Comments  */
    for (i=0; i<lenline; ++i) {
        if (line[i]=='!') {
            cstr[i] = '\0';
            break;
        }
        cstr[i] = line[i];
    }

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
    for (i=0; i<lenkname*6; i++) {
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

    /* AR  */
    kname[ 3*lenkname + 0 ] = 'A';
    kname[ 3*lenkname + 1 ] = 'R';
    kname[ 3*lenkname + 2 ] = ' ';

    /* HE  */
    kname[ 4*lenkname + 0 ] = 'H';
    kname[ 4*lenkname + 1 ] = 'E';
    kname[ 4*lenkname + 2 ] = ' ';

    /* C  */
    kname[ 5*lenkname + 0 ] = 'C';
    kname[ 5*lenkname + 1 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*11; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H2  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = ' ';

    /* O  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = 'H';
    kname[ 3*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 4*lenkname + 0 ] = 'H';
    kname[ 4*lenkname + 1 ] = '2';
    kname[ 4*lenkname + 2 ] = 'O';
    kname[ 4*lenkname + 3 ] = ' ';

    /* O2  */
    kname[ 5*lenkname + 0 ] = 'O';
    kname[ 5*lenkname + 1 ] = '2';
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

    /* AR  */
    kname[ 9*lenkname + 0 ] = 'A';
    kname[ 9*lenkname + 1 ] = 'R';
    kname[ 9*lenkname + 2 ] = ' ';

    /* HE  */
    kname[ 10*lenkname + 0 ] = 'H';
    kname[ 10*lenkname + 1 ] = 'E';
    kname[ 10*lenkname + 2 ] = ' ';

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
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
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

    for (int n=0; n<11; n++) {
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
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */
    W += c[9]*39.948000; /*AR */
    W += c[10]*4.002600; /*HE */

    for (id = 0; id < 11; ++id) {
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
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double YOW = 0;
    double tmp[11];

    for (int i = 0; i < 11; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 11; i++)
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
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */
    W += c[9]*39.948000; /*AR */
    W += c[10]*4.002600; /*HE */

    for (id = 0; id < 11; ++id) {
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
    double tmp[11];

    for (int i = 0; i < 11; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 11; i++)
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
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
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
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */
    W += c[9]*39.948000; /*AR */
    W += c[10]*4.002600; /*HE */

    for (id = 0; id < 11; ++id) {
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
    double tmp[11];

    for (int i = 0; i < 11; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 11; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 11; i++)
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

    for (int n=0; n<11; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<11; n++) {
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
    for (int i = 0; i < 11; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 11; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 11; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int i = 0; i < 11; i++)
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
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*1.007970*XWinv; 
    y[1] = x[1]*2.015940*XWinv; 
    y[2] = x[2]*15.999400*XWinv; 
    y[3] = x[3]*17.007370*XWinv; 
    y[4] = x[4]*18.015340*XWinv; 
    y[5] = x[5]*31.998800*XWinv; 
    y[6] = x[6]*33.006770*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*28.013400*XWinv; 
    y[9] = x[9]*39.948000*XWinv; 
    y[10] = x[10]*4.002600*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 11; ++id) {
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
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 11; ++id) {
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
    for (id = 0; id < 11; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 11; ++id) {
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
    CW += c[1]*2.015940; /*H2 */
    CW += c[2]*15.999400; /*O */
    CW += c[3]*17.007370; /*OH */
    CW += c[4]*18.015340; /*H2O */
    CW += c[5]*31.998800; /*O2 */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.013400; /*N2 */
    CW += c[9]*39.948000; /*AR */
    CW += c[10]*4.002600; /*HE */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*1.007970*CWinv; 
    y[1] = c[1]*2.015940*CWinv; 
    y[2] = c[2]*15.999400*CWinv; 
    y[3] = c[3]*17.007370*CWinv; 
    y[4] = c[4]*18.015340*CWinv; 
    y[5] = c[5]*31.998800*CWinv; 
    y[6] = c[6]*33.006770*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*28.013400*CWinv; 
    y[9] = c[9]*39.948000*CWinv; 
    y[10] = c[10]*4.002600*CWinv; 

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
    for (id = 0; id < 11; ++id) {
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
    for (id = 0; id < 11; ++id) {
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
    for (id = 0; id < 11; ++id) {
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
    for (id = 0; id < 11; ++id) {
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
    for (id = 0; id < 11; ++id) {
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
    for (id = 0; id < 11; ++id) {
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
    for (id = 0; id < 11; ++id) {
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
    cvms[1] *= 4.124383662212169e+07; /*H2 */
    cvms[2] *= 5.196763628636074e+06; /*O */
    cvms[3] *= 4.888768810227566e+06; /*OH */
    cvms[4] *= 4.615239012974499e+06; /*H2O */
    cvms[5] *= 2.598381814318037e+06; /*O2 */
    cvms[6] *= 2.519031701678171e+06; /*HO2 */
    cvms[7] *= 2.444384405113783e+06; /*H2O2 */
    cvms[8] *= 2.968047434442088e+06; /*N2 */
    cvms[9] *= 2.081333233203164e+06; /*AR */
    cvms[10] *= 2.077277269774646e+07; /*HE */
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
    cpms[1] *= 4.124383662212169e+07; /*H2 */
    cpms[2] *= 5.196763628636074e+06; /*O */
    cpms[3] *= 4.888768810227566e+06; /*OH */
    cpms[4] *= 4.615239012974499e+06; /*H2O */
    cpms[5] *= 2.598381814318037e+06; /*O2 */
    cpms[6] *= 2.519031701678171e+06; /*HO2 */
    cpms[7] *= 2.444384405113783e+06; /*H2O2 */
    cpms[8] *= 2.968047434442088e+06; /*N2 */
    cpms[9] *= 2.081333233203164e+06; /*AR */
    cpms[10] *= 2.077277269774646e+07; /*HE */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 11; i++)
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
    for (int i = 0; i < 11; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tc[5], h[11];

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
    }

    for (int n=0; n<11; n++) {
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
    for (int i = 0; i < 11; i++)
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
    for (int i = 0; i < 11; i++)
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
    sms[1] *= 4.124383662212169e+07; /*H2 */
    sms[2] *= 5.196763628636074e+06; /*O */
    sms[3] *= 4.888768810227566e+06; /*OH */
    sms[4] *= 4.615239012974499e+06; /*H2O */
    sms[5] *= 2.598381814318037e+06; /*O2 */
    sms[6] *= 2.519031701678171e+06; /*HO2 */
    sms[7] *= 2.444384405113783e+06; /*H2O2 */
    sms[8] *= 2.968047434442088e+06; /*N2 */
    sms[9] *= 2.081333233203164e+06; /*AR */
    sms[10] *= 2.077277269774646e+07; /*HE */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[11]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 11; ++id) {
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
    double cpor[11], tresult[11]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 11; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 11; i++)
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
    double cvor[11]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 11; ++id) {
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
    double cvor[11]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H */
    result += cvor[1]*y[1]*imw[1]; /*H2 */
    result += cvor[2]*y[2]*imw[2]; /*O */
    result += cvor[3]*y[3]*imw[3]; /*OH */
    result += cvor[4]*y[4]*imw[4]; /*H2O */
    result += cvor[5]*y[5]*imw[5]; /*O2 */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*N2 */
    result += cvor[9]*y[9]*imw[9]; /*AR */
    result += cvor[10]*y[10]*imw[10]; /*HE */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[11]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 11; ++id) {
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
    double hml[11], tmp[11]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 11; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 11; ++id) {
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
    double uml[11]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 11; ++id) {
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
    double ums[11]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H */
    result += y[1]*ums[1]*imw[1]; /*H2 */
    result += y[2]*ums[2]*imw[2]; /*O */
    result += y[3]*ums[3]*imw[3]; /*OH */
    result += y[4]*ums[4]*imw[4]; /*H2O */
    result += y[5]*ums[5]*imw[5]; /*O2 */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*N2 */
    result += y[9]*ums[9]*imw[9]; /*AR */
    result += y[10]*ums[10]*imw[10]; /*HE */

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
    double sor[11]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 11; ++id) {
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
    double sor[11]; /* temporary storage */
    double x[11]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(31.998800*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    x[9] = y[9]/(39.948000*YOW); 
    x[10] = y[10]/(4.002600*YOW); 
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
    double gort[11]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 11; ++id) {
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
    double gort[11]; /* temporary storage */
    double x[11]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(31.998800*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    x[9] = y[9]/(39.948000*YOW); 
    x[10] = y[10]/(4.002600*YOW); 
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
    double aort[11]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 11; ++id) {
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
    double aort[11]; /* temporary storage */
    double x[11]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(31.998800*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    x[9] = y[9]/(39.948000*YOW); 
    x[10] = y[10]/(4.002600*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 11; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 11; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[11]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
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

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 11; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[11]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 11; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 11; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[11]; /*temporary storage */
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

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 11; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[11*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<11; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<11*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[11]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 11; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 11; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 11; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 11; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 29; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)
{
    int id; /*loop counter */
    double c[11]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 11; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[11]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
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

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[11]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 11; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[11]; /*temporary storage */
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

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[11]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 11; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
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
    for (id = 0; id < 11 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 0 * kd + 0 ] += -1 ;
    nuki[ 5 * kd + 0 ] += -1 ;
    nuki[ 6 * kd + 0 ] += +1 ;

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    nuki[ 7 * kd + 1 ] += -1 ;
    nuki[ 3 * kd + 1 ] += +1 ;
    nuki[ 3 * kd + 1 ] += +1 ;

    /*reaction 3: H2 + M <=> H + H + M */
    nuki[ 1 * kd + 2 ] += -1 ;
    nuki[ 0 * kd + 2 ] += +1 ;
    nuki[ 0 * kd + 2 ] += +1 ;

    /*reaction 4: O + O + M <=> O2 + M */
    nuki[ 2 * kd + 3 ] += -1 ;
    nuki[ 2 * kd + 3 ] += -1 ;
    nuki[ 5 * kd + 3 ] += +1 ;

    /*reaction 5: O + H + M <=> OH + M */
    nuki[ 2 * kd + 4 ] += -1 ;
    nuki[ 0 * kd + 4 ] += -1 ;
    nuki[ 3 * kd + 4 ] += +1 ;

    /*reaction 6: H2O + M <=> H + OH + M */
    nuki[ 4 * kd + 5 ] += -1 ;
    nuki[ 0 * kd + 5 ] += +1 ;
    nuki[ 3 * kd + 5 ] += +1 ;

    /*reaction 7: O + OH + M <=> HO2 + M */
    nuki[ 2 * kd + 6 ] += -1 ;
    nuki[ 3 * kd + 6 ] += -1 ;
    nuki[ 6 * kd + 6 ] += +1 ;

    /*reaction 8: H + O2 <=> O + OH */
    nuki[ 0 * kd + 7 ] += -1 ;
    nuki[ 5 * kd + 7 ] += -1 ;
    nuki[ 2 * kd + 7 ] += +1 ;
    nuki[ 3 * kd + 7 ] += +1 ;

    /*reaction 9: O + H2 <=> H + OH */
    nuki[ 2 * kd + 8 ] += -1 ;
    nuki[ 1 * kd + 8 ] += -1 ;
    nuki[ 0 * kd + 8 ] += +1 ;
    nuki[ 3 * kd + 8 ] += +1 ;

    /*reaction 10: O + H2 <=> H + OH */
    nuki[ 2 * kd + 9 ] += -1 ;
    nuki[ 1 * kd + 9 ] += -1 ;
    nuki[ 0 * kd + 9 ] += +1 ;
    nuki[ 3 * kd + 9 ] += +1 ;

    /*reaction 11: H2 + OH <=> H2O + H */
    nuki[ 1 * kd + 10 ] += -1 ;
    nuki[ 3 * kd + 10 ] += -1 ;
    nuki[ 4 * kd + 10 ] += +1 ;
    nuki[ 0 * kd + 10 ] += +1 ;

    /*reaction 12: OH + OH <=> O + H2O */
    nuki[ 3 * kd + 11 ] += -1 ;
    nuki[ 3 * kd + 11 ] += -1 ;
    nuki[ 2 * kd + 11 ] += +1 ;
    nuki[ 4 * kd + 11 ] += +1 ;

    /*reaction 13: H2 + AR <=> H + H + AR */
    nuki[ 1 * kd + 12 ] += -1 ;
    nuki[ 9 * kd + 12 ] += -1 ;
    nuki[ 0 * kd + 12 ] += +1 ;
    nuki[ 0 * kd + 12 ] += +1 ;
    nuki[ 9 * kd + 12 ] += +1 ;

    /*reaction 14: H2 + HE <=> H + H + HE */
    nuki[ 1 * kd + 13 ] += -1 ;
    nuki[ 10 * kd + 13 ] += -1 ;
    nuki[ 0 * kd + 13 ] += +1 ;
    nuki[ 0 * kd + 13 ] += +1 ;
    nuki[ 10 * kd + 13 ] += +1 ;

    /*reaction 15: O + O + AR <=> O2 + AR */
    nuki[ 2 * kd + 14 ] += -1 ;
    nuki[ 2 * kd + 14 ] += -1 ;
    nuki[ 9 * kd + 14 ] += -1 ;
    nuki[ 5 * kd + 14 ] += +1 ;
    nuki[ 9 * kd + 14 ] += +1 ;

    /*reaction 16: O + O + HE <=> O2 + HE */
    nuki[ 2 * kd + 15 ] += -1 ;
    nuki[ 2 * kd + 15 ] += -1 ;
    nuki[ 10 * kd + 15 ] += -1 ;
    nuki[ 5 * kd + 15 ] += +1 ;
    nuki[ 10 * kd + 15 ] += +1 ;

    /*reaction 17: H2O + H2O <=> H + OH + H2O */
    nuki[ 4 * kd + 16 ] += -1 ;
    nuki[ 4 * kd + 16 ] += -1 ;
    nuki[ 0 * kd + 16 ] += +1 ;
    nuki[ 3 * kd + 16 ] += +1 ;
    nuki[ 4 * kd + 16 ] += +1 ;

    /*reaction 18: HO2 + H <=> H2 + O2 */
    nuki[ 6 * kd + 17 ] += -1 ;
    nuki[ 0 * kd + 17 ] += -1 ;
    nuki[ 1 * kd + 17 ] += +1 ;
    nuki[ 5 * kd + 17 ] += +1 ;

    /*reaction 19: HO2 + H <=> OH + OH */
    nuki[ 6 * kd + 18 ] += -1 ;
    nuki[ 0 * kd + 18 ] += -1 ;
    nuki[ 3 * kd + 18 ] += +1 ;
    nuki[ 3 * kd + 18 ] += +1 ;

    /*reaction 20: HO2 + O <=> O2 + OH */
    nuki[ 6 * kd + 19 ] += -1 ;
    nuki[ 2 * kd + 19 ] += -1 ;
    nuki[ 5 * kd + 19 ] += +1 ;
    nuki[ 3 * kd + 19 ] += +1 ;

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    nuki[ 6 * kd + 20 ] += -1 ;
    nuki[ 3 * kd + 20 ] += -1 ;
    nuki[ 4 * kd + 20 ] += +1 ;
    nuki[ 5 * kd + 20 ] += +1 ;

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 6 * kd + 21 ] += -1 ;
    nuki[ 6 * kd + 21 ] += -1 ;
    nuki[ 7 * kd + 21 ] += +1 ;
    nuki[ 5 * kd + 21 ] += +1 ;

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 6 * kd + 22 ] += -1 ;
    nuki[ 6 * kd + 22 ] += -1 ;
    nuki[ 7 * kd + 22 ] += +1 ;
    nuki[ 5 * kd + 22 ] += +1 ;

    /*reaction 24: H2O2 + H <=> H2O + OH */
    nuki[ 7 * kd + 23 ] += -1 ;
    nuki[ 0 * kd + 23 ] += -1 ;
    nuki[ 4 * kd + 23 ] += +1 ;
    nuki[ 3 * kd + 23 ] += +1 ;

    /*reaction 25: H2O2 + H <=> HO2 + H2 */
    nuki[ 7 * kd + 24 ] += -1 ;
    nuki[ 0 * kd + 24 ] += -1 ;
    nuki[ 6 * kd + 24 ] += +1 ;
    nuki[ 1 * kd + 24 ] += +1 ;

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    nuki[ 7 * kd + 25 ] += -1 ;
    nuki[ 2 * kd + 25 ] += -1 ;
    nuki[ 3 * kd + 25 ] += +1 ;
    nuki[ 6 * kd + 25 ] += +1 ;

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    nuki[ 7 * kd + 26 ] += -1 ;
    nuki[ 3 * kd + 26 ] += -1 ;
    nuki[ 6 * kd + 26 ] += +1 ;
    nuki[ 4 * kd + 26 ] += +1 ;

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    nuki[ 7 * kd + 27 ] += -1 ;
    nuki[ 3 * kd + 27 ] += -1 ;
    nuki[ 6 * kd + 27 ] += +1 ;
    nuki[ 4 * kd + 27 ] += +1 ;

    /*reaction 29: HO2 + H <=> O + H2O */
    nuki[ 6 * kd + 28 ] += -1 ;
    nuki[ 0 * kd + 28 ] += -1 ;
    nuki[ 2 * kd + 28 ] += +1 ;
    nuki[ 4 * kd + 28 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 11; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 0 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 0 ] = 2; /*H */

    /*O */
    ncf[ 2 * kd + 1 ] = 1; /*O */

    /*OH */
    ncf[ 3 * kd + 1 ] = 1; /*O */
    ncf[ 3 * kd + 0 ] = 1; /*H */

    /*H2O */
    ncf[ 4 * kd + 0 ] = 2; /*H */
    ncf[ 4 * kd + 1 ] = 1; /*O */

    /*O2 */
    ncf[ 5 * kd + 1 ] = 2; /*O */

    /*HO2 */
    ncf[ 6 * kd + 0 ] = 1; /*H */
    ncf[ 6 * kd + 1 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 0 ] = 2; /*H */
    ncf[ 7 * kd + 1 ] = 2; /*O */

    /*N2 */
    ncf[ 8 * kd + 2 ] = 2; /*N */

    /*AR */
    ncf[ 9 * kd + 3 ] = 1; /*AR */

    /*HE */
    ncf[ 10 * kd + 4 ] = 1; /*HE */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)
{
    for (int i=0; i<29; ++i) {
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
    double gort[11]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2O + M <=> H + OH + M */
    eqcon[5] *= 1e-06; 

    /*reaction 7: O + OH + M <=> HO2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 <=> O + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2 <=> H + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H2 + AR <=> H + H + AR */
    eqcon[12] *= 1e-06; 

    /*reaction 14: H2 + HE <=> H + H + HE */
    eqcon[13] *= 1e-06; 

    /*reaction 15: O + O + AR <=> O2 + AR */
    eqcon[14] *= 1e+06; 

    /*reaction 16: O + O + HE <=> O2 + HE */
    eqcon[15] *= 1e+06; 

    /*reaction 17: H2O + H2O <=> H + OH + H2O */
    eqcon[16] *= 1e-06; 

    /*reaction 18: HO2 + H <=> H2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + O <=> O2 + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> H2O + OH */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> HO2 + H2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HO2 + H <=> O + H2O */
    /*eqcon[28] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[11]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2O + M <=> H + OH + M */
    eqcon[5] *= 1e-06; 

    /*reaction 7: O + OH + M <=> HO2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 <=> O + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2 <=> H + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H2 + AR <=> H + H + AR */
    eqcon[12] *= 1e-06; 

    /*reaction 14: H2 + HE <=> H + H + HE */
    eqcon[13] *= 1e-06; 

    /*reaction 15: O + O + AR <=> O2 + AR */
    eqcon[14] *= 1e+06; 

    /*reaction 16: O + O + HE <=> O2 + HE */
    eqcon[15] *= 1e+06; 

    /*reaction 17: H2O + H2O <=> H + OH + H2O */
    eqcon[16] *= 1e-06; 

    /*reaction 18: HO2 + H <=> H2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + O <=> O2 + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> H2O + OH */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> HO2 + H2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HO2 + H <=> O + H2O */
    /*eqcon[28] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[11]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2O + M <=> H + OH + M */
    eqcon[5] *= 1e-06; 

    /*reaction 7: O + OH + M <=> HO2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 <=> O + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2 <=> H + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H2 + AR <=> H + H + AR */
    eqcon[12] *= 1e-06; 

    /*reaction 14: H2 + HE <=> H + H + HE */
    eqcon[13] *= 1e-06; 

    /*reaction 15: O + O + AR <=> O2 + AR */
    eqcon[14] *= 1e+06; 

    /*reaction 16: O + O + HE <=> O2 + HE */
    eqcon[15] *= 1e+06; 

    /*reaction 17: H2O + H2O <=> H + OH + H2O */
    eqcon[16] *= 1e-06; 

    /*reaction 18: HO2 + H <=> H2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + O <=> O2 + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> H2O + OH */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> HO2 + H2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HO2 + H <=> O + H2O */
    /*eqcon[28] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[11]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2O + M <=> H + OH + M */
    eqcon[5] *= 1e-06; 

    /*reaction 7: O + OH + M <=> HO2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 <=> O + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2 <=> H + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H2 + AR <=> H + H + AR */
    eqcon[12] *= 1e-06; 

    /*reaction 14: H2 + HE <=> H + H + HE */
    eqcon[13] *= 1e-06; 

    /*reaction 15: O + O + AR <=> O2 + AR */
    eqcon[14] *= 1e+06; 

    /*reaction 16: O + O + HE <=> O2 + HE */
    eqcon[15] *= 1e+06; 

    /*reaction 17: H2O + H2O <=> H + OH + H2O */
    eqcon[16] *= 1e-06; 

    /*reaction 18: HO2 + H <=> H2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + O <=> O2 + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> H2O + OH */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> HO2 + H2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HO2 + H <=> O + H2O */
    /*eqcon[28] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[11]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2O + M <=> H + OH + M */
    eqcon[5] *= 1e-06; 

    /*reaction 7: O + OH + M <=> HO2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 <=> O + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2 <=> H + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H2 + AR <=> H + H + AR */
    eqcon[12] *= 1e-06; 

    /*reaction 14: H2 + HE <=> H + H + HE */
    eqcon[13] *= 1e-06; 

    /*reaction 15: O + O + AR <=> O2 + AR */
    eqcon[14] *= 1e+06; 

    /*reaction 16: O + O + HE <=> O2 + HE */
    eqcon[15] *= 1e+06; 

    /*reaction 17: H2O + H2O <=> H + OH + H2O */
    eqcon[16] *= 1e-06; 

    /*reaction 18: HO2 + H <=> H2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + O <=> O2 + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> H2O + OH */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> HO2 + H2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HO2 + H <=> O + H2O */
    /*eqcon[28] *= 1;  */
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[29];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[29];
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

    double qdot, q_f[29], q_r[29];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 11; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[0] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[2]-q_r[2];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;

    qdot = q_f[3]-q_r[3];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[0] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[10] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[10] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[0] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[4] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[21]-q_r[21];
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[23]-q_r[23];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    return;
}

void comp_k_f(double * restrict tc, double invT, double * restrict k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<29; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double * restrict tc, double invT, double * restrict Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[11];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[0] + g_RT[5] - g_RT[6];
    Kc[1] = -g_RT[3] - g_RT[3] + g_RT[7];
    Kc[2] = -g_RT[0] - g_RT[0] + g_RT[1];
    Kc[3] = g_RT[2] + g_RT[2] - g_RT[5];
    Kc[4] = g_RT[0] + g_RT[2] - g_RT[3];
    Kc[5] = -g_RT[0] - g_RT[3] + g_RT[4];
    Kc[6] = g_RT[2] + g_RT[3] - g_RT[6];
    Kc[7] = g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5];
    Kc[8] = -g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3];
    Kc[9] = -g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3];
    Kc[10] = -g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4];
    Kc[11] = -g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4];
    Kc[12] = -g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9];
    Kc[13] = -g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10];
    Kc[14] = g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9];
    Kc[15] = g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10];
    Kc[16] = -g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4];
    Kc[17] = g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6];
    Kc[18] = g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6];
    Kc[19] = g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6];
    Kc[20] = g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6];
    Kc[21] = -g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7];
    Kc[22] = -g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7];
    Kc[23] = g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7];
    Kc[24] = g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7];
    Kc[25] = g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7];
    Kc[26] = g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7];
    Kc[27] = g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7];
    Kc[28] = g_RT[0] - g_RT[2] - g_RT[4] + g_RT[6];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<29; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refC;
    Kc[2] *= refC;
    Kc[3] *= refCinv;
    Kc[4] *= refCinv;
    Kc[5] *= refC;
    Kc[6] *= refCinv;
    Kc[12] *= refC;
    Kc[13] *= refC;
    Kc[14] *= refCinv;
    Kc[15] *= refCinv;
    Kc[16] *= refC;

    return;
}

void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[0]*sc[5];
    qr[0] = sc[6];

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    qf[1] = sc[7];
    qr[1] = sc[3]*sc[3];

    /*reaction 3: H2 + M <=> H + H + M */
    qf[2] = sc[1];
    qr[2] = sc[0]*sc[0];

    /*reaction 4: O + O + M <=> O2 + M */
    qf[3] = sc[2]*sc[2];
    qr[3] = sc[5];

    /*reaction 5: O + H + M <=> OH + M */
    qf[4] = sc[0]*sc[2];
    qr[4] = sc[3];

    /*reaction 6: H2O + M <=> H + OH + M */
    qf[5] = sc[4];
    qr[5] = sc[0]*sc[3];

    /*reaction 7: O + OH + M <=> HO2 + M */
    qf[6] = sc[2]*sc[3];
    qr[6] = sc[6];

    /*reaction 8: H + O2 <=> O + OH */
    qf[7] = sc[0]*sc[5];
    qr[7] = sc[2]*sc[3];

    /*reaction 9: O + H2 <=> H + OH */
    qf[8] = sc[1]*sc[2];
    qr[8] = sc[0]*sc[3];

    /*reaction 10: O + H2 <=> H + OH */
    qf[9] = sc[1]*sc[2];
    qr[9] = sc[0]*sc[3];

    /*reaction 11: H2 + OH <=> H2O + H */
    qf[10] = sc[1]*sc[3];
    qr[10] = sc[0]*sc[4];

    /*reaction 12: OH + OH <=> O + H2O */
    qf[11] = sc[3]*sc[3];
    qr[11] = sc[2]*sc[4];

    /*reaction 13: H2 + AR <=> H + H + AR */
    qf[12] = sc[1]*sc[9];
    qr[12] = sc[0]*sc[0]*sc[9];

    /*reaction 14: H2 + HE <=> H + H + HE */
    qf[13] = sc[1]*sc[10];
    qr[13] = sc[0]*sc[0]*sc[10];

    /*reaction 15: O + O + AR <=> O2 + AR */
    qf[14] = sc[2]*sc[2]*sc[9];
    qr[14] = sc[5]*sc[9];

    /*reaction 16: O + O + HE <=> O2 + HE */
    qf[15] = sc[2]*sc[2]*sc[10];
    qr[15] = sc[5]*sc[10];

    /*reaction 17: H2O + H2O <=> H + OH + H2O */
    qf[16] = sc[4]*sc[4];
    qr[16] = sc[0]*sc[3]*sc[4];

    /*reaction 18: HO2 + H <=> H2 + O2 */
    qf[17] = sc[0]*sc[6];
    qr[17] = sc[1]*sc[5];

    /*reaction 19: HO2 + H <=> OH + OH */
    qf[18] = sc[0]*sc[6];
    qr[18] = sc[3]*sc[3];

    /*reaction 20: HO2 + O <=> O2 + OH */
    qf[19] = sc[2]*sc[6];
    qr[19] = sc[3]*sc[5];

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    qf[20] = sc[3]*sc[6];
    qr[20] = sc[4]*sc[5];

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    qf[21] = sc[6]*sc[6];
    qr[21] = sc[5]*sc[7];

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    qf[22] = sc[6]*sc[6];
    qr[22] = sc[5]*sc[7];

    /*reaction 24: H2O2 + H <=> H2O + OH */
    qf[23] = sc[0]*sc[7];
    qr[23] = sc[3]*sc[4];

    /*reaction 25: H2O2 + H <=> HO2 + H2 */
    qf[24] = sc[0]*sc[7];
    qr[24] = sc[1]*sc[6];

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    qf[25] = sc[2]*sc[7];
    qr[25] = sc[3]*sc[6];

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    qf[26] = sc[3]*sc[7];
    qr[26] = sc[4]*sc[6];

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    qf[27] = sc[3]*sc[7];
    qr[27] = sc[4]*sc[6];

    /*reaction 29: HO2 + H <=> O + H2O */
    qf[28] = sc[0]*sc[6];
    qr[28] = sc[2]*sc[4];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 11; ++i) {
        mixture += sc[i];
    }

    double Corr[29];
    for (int i = 0; i < 29; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[2];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[1] + (TB[0][1] - 1)*sc[4] + (TB[0][2] - 1)*sc[5] + (TB[0][3] - 1)*sc[9] + (TB[0][4] - 1)*sc[10];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[4] + (TB[1][1] - 1)*sc[8] + (TB[1][2] - 1)*sc[5] + (TB[1][3] - 1)*sc[10] + (TB[1][4] - 1)*sc[7] + (TB[1][5] - 1)*sc[1];
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

    /* simple three-body correction */
    {
        double alpha;
        alpha = mixture + (TB[2][0] - 1)*sc[1] + (TB[2][1] - 1)*sc[4] + (TB[2][2] - 1)*sc[9] + (TB[2][3] - 1)*sc[10];
        Corr[2] = alpha;
        alpha = mixture + (TB[3][0] - 1)*sc[1] + (TB[3][1] - 1)*sc[4] + (TB[3][2] - 1)*sc[9] + (TB[3][3] - 1)*sc[10];
        Corr[3] = alpha;
        alpha = mixture + (TB[4][0] - 1)*sc[1] + (TB[4][1] - 1)*sc[4] + (TB[4][2] - 1)*sc[9] + (TB[4][3] - 1)*sc[10];
        Corr[4] = alpha;
        alpha = mixture + (TB[5][0] - 1)*sc[1] + (TB[5][1] - 1)*sc[4] + (TB[5][2] - 1)*sc[10] + (TB[5][3] - 1)*sc[8] + (TB[5][4] - 1)*sc[5];
        Corr[5] = alpha;
        alpha = mixture + (TB[6][0] - 1)*sc[1] + (TB[6][1] - 1)*sc[4] + (TB[6][2] - 1)*sc[9] + (TB[6][3] - 1)*sc[10];
        Corr[6] = alpha;
    }

    for (int i=0; i<29; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[29*npt], Kc_s[29*npt], mixture[npt], g_RT[11*npt];
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

    for (int n=0; n<11; n++) {
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
    }
}

void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[11];
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

        Kc_s[0*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[1*npt+i] = refC * exp((g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[2*npt+i] = refC * exp((g_RT[1*npt+i]) - (g_RT[0*npt+i] + g_RT[0*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[2*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[3*npt+i]));
        Kc_s[5*npt+i] = refC * exp((g_RT[4*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[7*npt+i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[9*npt+i] = exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[10*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[0*npt+i] + g_RT[4*npt+i]));
        Kc_s[11*npt+i] = exp((g_RT[3*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
        Kc_s[12*npt+i] = refC * exp((g_RT[1*npt+i] + g_RT[9*npt+i]) - (g_RT[0*npt+i] + g_RT[0*npt+i] + g_RT[9*npt+i]));
        Kc_s[13*npt+i] = refC * exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[0*npt+i] + g_RT[0*npt+i] + g_RT[10*npt+i]));
        Kc_s[14*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[2*npt+i] + g_RT[9*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[15*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[2*npt+i] + g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[10*npt+i]));
        Kc_s[16*npt+i] = refC * exp((g_RT[4*npt+i] + g_RT[4*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[0*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[0*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[2*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[3*npt+i] + g_RT[6*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[0*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[0*npt+i] + g_RT[7*npt+i]) - (g_RT[1*npt+i] + g_RT[6*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[0*npt+i] + g_RT[6*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
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
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[1*npt+i] + (TB[0][1] - 1)*sc[4*npt+i] + (TB[0][2] - 1)*sc[5*npt+i] + (TB[0][3] - 1)*sc[9*npt+i] + (TB[0][4] - 1)*sc[10*npt+i];
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
        phi_r = sc[6*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[4*npt+i] + (TB[1][1] - 1)*sc[8*npt+i] + (TB[1][2] - 1)*sc[5*npt+i] + (TB[1][3] - 1)*sc[10*npt+i] + (TB[1][4] - 1)*sc[7*npt+i] + (TB[1][5] - 1)*sc[1*npt+i];
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
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 3: H2 + M <=> H + H + M */
        phi_f = sc[1*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[1*npt+i] + (TB[2][1] - 1)*sc[4*npt+i] + (TB[2][2] - 1)*sc[9*npt+i] + (TB[2][3] - 1)*sc[10*npt+i];
        k_f = alpha * k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;

        /*reaction 4: O + O + M <=> O2 + M */
        phi_f = sc[2*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[1*npt+i] + (TB[3][1] - 1)*sc[4*npt+i] + (TB[3][2] - 1)*sc[9*npt+i] + (TB[3][3] - 1)*sc[10*npt+i];
        k_f = alpha * k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 5: O + H + M <=> OH + M */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[1*npt+i] + (TB[4][1] - 1)*sc[4*npt+i] + (TB[4][2] - 1)*sc[9*npt+i] + (TB[4][3] - 1)*sc[10*npt+i];
        k_f = alpha * k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 6: H2O + M <=> H + OH + M */
        phi_f = sc[4*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[1*npt+i] + (TB[5][1] - 1)*sc[4*npt+i] + (TB[5][2] - 1)*sc[10*npt+i] + (TB[5][3] - 1)*sc[8*npt+i] + (TB[5][4] - 1)*sc[5*npt+i];
        k_f = alpha * k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;

        /*reaction 7: O + OH + M <=> HO2 + M */
        phi_f = sc[2*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[1*npt+i] + (TB[6][1] - 1)*sc[4*npt+i] + (TB[6][2] - 1)*sc[9*npt+i] + (TB[6][3] - 1)*sc[10*npt+i];
        k_f = alpha * k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 8: H + O2 <=> O + OH */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[3*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 9: O + H2 <=> H + OH */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 10: O + H2 <=> H + OH */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 11: H2 + OH <=> H2O + H */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[4*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 12: OH + OH <=> O + H2O */
        phi_f = sc[3*npt+i]*sc[3*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 13: H2 + AR <=> H + H + AR */
        phi_f = sc[1*npt+i]*sc[9*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i]*sc[9*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 14: H2 + HE <=> H + H + HE */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i]*sc[10*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 15: O + O + AR <=> O2 + AR */
        phi_f = sc[2*npt+i]*sc[2*npt+i]*sc[9*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 16: O + O + HE <=> O2 + HE */
        phi_f = sc[2*npt+i]*sc[2*npt+i]*sc[10*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[10*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 17: H2O + H2O <=> H + OH + H2O */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 18: HO2 + H <=> H2 + O2 */
        phi_f = sc[0*npt+i]*sc[6*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 19: HO2 + H <=> OH + OH */
        phi_f = sc[0*npt+i]*sc[6*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 20: HO2 + O <=> O2 + OH */
        phi_f = sc[2*npt+i]*sc[6*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 21: HO2 + OH <=> H2O + O2 */
        phi_f = sc[3*npt+i]*sc[6*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[7*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[7*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 24: H2O2 + H <=> H2O + OH */
        phi_f = sc[0*npt+i]*sc[7*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 25: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[0*npt+i]*sc[7*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[6*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 26: H2O2 + O <=> OH + HO2 */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 27: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 28: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 29: HO2 + H <=> O + H2O */
        phi_f = sc[0*npt+i]*sc[6*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
    }
}

/*compute the reaction Jacobian */
void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)
{
    double c[11];

    for (int k=0; k<11; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    for (int k=0; k<11; k++) {
        J[132+k] *= 1.e-6;
    }

    /* dTdot/d[X] */
    for (int k=0; k<11; k++) {
        J[k*12+11] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<144; i++) {
        J[i] = 0.0;
    }

    double wdot[11];
    for (int k=0; k<11; k++) {
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
    for (int k = 0; k < 11; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[11];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[11];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[11];
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
    alpha = mixture + (TB[0][0] - 1)*sc[1] + (TB[0][1] - 1)*sc[4] + (TB[0][2] - 1)*sc[5] + (TB[0][3] - 1)*sc[9] + (TB[0][4] - 1)*sc[10];
    /* forward */
    phi_f = sc[0]*sc[5];
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
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[5] -= dqdci;                /* dwdot[O2]/d[H] */
        J[6] += dqdci;                /* dwdot[HO2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[12] -= dqdci;               /* dwdot[H]/d[H2] */
        J[17] -= dqdci;               /* dwdot[O2]/d[H2] */
        J[18] += dqdci;               /* dwdot[HO2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[48] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[53] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[54] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (TB[0][2] - 1)*dcdc_fac + k_f*sc[0];
        J[60] -= dqdci;               /* dwdot[H]/d[O2] */
        J[65] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[66] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[72] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[77] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[78] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[AR] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[108] -= dqdci;              /* dwdot[H]/d[AR] */
        J[113] -= dqdci;              /* dwdot[O2]/d[AR] */
        J[114] += dqdci;              /* dwdot[HO2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[0][4] - 1)*dcdc_fac;
        J[120] -= dqdci;              /* dwdot[H]/d[HE] */
        J[125] -= dqdci;              /* dwdot[O2]/d[HE] */
        J[126] += dqdci;              /* dwdot[HO2]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[5];
        dqdc[1] = TB[0][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[0][1]*dcdc_fac;
        dqdc[5] = TB[0][2]*dcdc_fac + k_f*sc[0];
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[0][3]*dcdc_fac;
        dqdc[10] = TB[0][4]*dcdc_fac;
        for (int k=0; k<11; k++) {
            J[12*k+0] -= dqdc[k];
            J[12*k+5] -= dqdc[k];
            J[12*k+6] += dqdc[k];
        }
    }
    J[132] -= dqdT; /* dwdot[H]/dT */
    J[137] -= dqdT; /* dwdot[O2]/dT */
    J[138] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[4] + (TB[1][1] - 1)*sc[8] + (TB[1][2] - 1)*sc[5] + (TB[1][3] - 1)*sc[10] + (TB[1][4] - 1)*sc[7] + (TB[1][5] - 1)*sc[1];
    /* forward */
    phi_f = sc[7];
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
    phi_r = sc[3]*sc[3];
    Kc = refC * exp(-g_RT[3] - g_RT[3] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2*h_RT[3]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[3] += 2 * q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[1][5] - 1)*dcdc_fac;
        J[15] += 2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[19] -= dqdci;               /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*2*sc[3];
        J[39] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
        J[43] -= dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[51] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[55] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (TB[1][2] - 1)*dcdc_fac;
        J[63] += 2 * dqdci;           /* dwdot[OH]/d[O2] */
        J[67] -= dqdci;               /* dwdot[H2O2]/d[O2] */
        /* d()/d[H2O2] */
        dqdci = (TB[1][4] - 1)*dcdc_fac + k_f;
        J[87] += 2 * dqdci;           /* dwdot[OH]/d[H2O2] */
        J[91] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[N2] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[99] += 2 * dqdci;           /* dwdot[OH]/d[N2] */
        J[103] -= dqdci;              /* dwdot[H2O2]/d[N2] */
        /* d()/d[HE] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[123] += 2 * dqdci;          /* dwdot[OH]/d[HE] */
        J[127] -= dqdci;              /* dwdot[H2O2]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[1][5]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac - k_r*2*sc[3];
        dqdc[4] = TB[1][0]*dcdc_fac;
        dqdc[5] = TB[1][2]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[1][4]*dcdc_fac + k_f;
        dqdc[8] = TB[1][1]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = TB[1][3]*dcdc_fac;
        for (int k=0; k<11; k++) {
            J[12*k+3] += 2 * dqdc[k];
            J[12*k+7] -= dqdc[k];
        }
    }
    J[135] += 2 * dqdT; /* dwdot[OH]/dT */
    J[139] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[1] + (TB[2][1] - 1)*sc[4] + (TB[2][2] - 1)*sc[9] + (TB[2][3] - 1)*sc[10];
    /* forward */
    phi_f = sc[1];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[0];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1]) + (2*h_RT[0]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*2*sc[0];
        J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
        J[1] -= dqdci;                /* dwdot[H2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*q_nocor + k_f;
        J[12] += 2 * dqdci;           /* dwdot[H]/d[H2] */
        J[13] -= dqdci;               /* dwdot[H2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*q_nocor;
        J[48] += 2 * dqdci;           /* dwdot[H]/d[H2O] */
        J[49] -= dqdci;               /* dwdot[H2]/d[H2O] */
        /* d()/d[AR] */
        dqdci = (TB[2][2] - 1)*q_nocor;
        J[108] += 2 * dqdci;          /* dwdot[H]/d[AR] */
        J[109] -= dqdci;              /* dwdot[H2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[2][3] - 1)*q_nocor;
        J[120] += 2 * dqdci;          /* dwdot[H]/d[HE] */
        J[121] -= dqdci;              /* dwdot[H2]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac - k_r*2*sc[0];
        dqdc[1] = TB[2][0] + k_f;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[2][1];
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[2][2];
        dqdc[10] = TB[2][3];
        for (int k=0; k<11; k++) {
            J[12*k+0] += 2 * dqdc[k];
            J[12*k+1] -= dqdc[k];
        }
    }
    J[132] += 2 * dqdT; /* dwdot[H]/dT */
    J[133] -= dqdT; /* dwdot[H2]/dT */

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[1] + (TB[3][1] - 1)*sc[4] + (TB[3][2] - 1)*sc[9] + (TB[3][3] - 1)*sc[10];
    /* forward */
    phi_f = sc[2]*sc[2];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[2]) + (h_RT[5]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*q_nocor;
        J[14] += -2 * dqdci;          /* dwdot[O]/d[H2] */
        J[17] += dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2*sc[2];
        J[26] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[29] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*q_nocor;
        J[50] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[62] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        J[65] += dqdci;               /* dwdot[O2]/d[O2] */
        /* d()/d[AR] */
        dqdci = (TB[3][2] - 1)*q_nocor;
        J[110] += -2 * dqdci;         /* dwdot[O]/d[AR] */
        J[113] += dqdci;              /* dwdot[O2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[3][3] - 1)*q_nocor;
        J[122] += -2 * dqdci;         /* dwdot[O]/d[HE] */
        J[125] += dqdci;              /* dwdot[O2]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[3][0];
        dqdc[2] = dcdc_fac + k_f*2*sc[2];
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[3][1];
        dqdc[5] = dcdc_fac - k_r;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[3][2];
        dqdc[10] = TB[3][3];
        for (int k=0; k<11; k++) {
            J[12*k+2] += -2 * dqdc[k];
            J[12*k+5] += dqdc[k];
        }
    }
    J[134] += -2 * dqdT; /* dwdot[O]/dT */
    J[137] += dqdT; /* dwdot[O2]/dT */

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[1] + (TB[4][1] - 1)*sc[4] + (TB[4][2] - 1)*sc[9] + (TB[4][3] - 1)*sc[10];
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[2] -= dqdci;                /* dwdot[O]/d[H] */
        J[3] += dqdci;                /* dwdot[OH]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*q_nocor;
        J[12] -= dqdci;               /* dwdot[H]/d[H2] */
        J[14] -= dqdci;               /* dwdot[O]/d[H2] */
        J[15] += dqdci;               /* dwdot[OH]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[0];
        J[24] -= dqdci;               /* dwdot[H]/d[O] */
        J[26] -= dqdci;               /* dwdot[O]/d[O] */
        J[27] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[36] -= dqdci;               /* dwdot[H]/d[OH] */
        J[38] -= dqdci;               /* dwdot[O]/d[OH] */
        J[39] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*q_nocor;
        J[48] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[50] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[51] += dqdci;               /* dwdot[OH]/d[H2O] */
        /* d()/d[AR] */
        dqdci = (TB[4][2] - 1)*q_nocor;
        J[108] -= dqdci;              /* dwdot[H]/d[AR] */
        J[110] -= dqdci;              /* dwdot[O]/d[AR] */
        J[111] += dqdci;              /* dwdot[OH]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[4][3] - 1)*q_nocor;
        J[120] -= dqdci;              /* dwdot[H]/d[HE] */
        J[122] -= dqdci;              /* dwdot[O]/d[HE] */
        J[123] += dqdci;              /* dwdot[OH]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[2];
        dqdc[1] = TB[4][0];
        dqdc[2] = dcdc_fac + k_f*sc[0];
        dqdc[3] = dcdc_fac - k_r;
        dqdc[4] = TB[4][1];
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[4][2];
        dqdc[10] = TB[4][3];
        for (int k=0; k<11; k++) {
            J[12*k+0] -= dqdc[k];
            J[12*k+2] -= dqdc[k];
            J[12*k+3] += dqdc[k];
        }
    }
    J[132] -= dqdT; /* dwdot[H]/dT */
    J[134] -= dqdT; /* dwdot[O]/dT */
    J[135] += dqdT; /* dwdot[OH]/dT */

    /*reaction 6: H2O + M <=> H + OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[1] + (TB[5][1] - 1)*sc[4] + (TB[5][2] - 1)*sc[10] + (TB[5][3] - 1)*sc[8] + (TB[5][4] - 1)*sc[5];
    /* forward */
    phi_f = sc[4];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4]) + (h_RT[0] + h_RT[3]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*sc[3];
        J[0] += dqdci;                /* dwdot[H]/d[H] */
        J[3] += dqdci;                /* dwdot[OH]/d[H] */
        J[4] -= dqdci;                /* dwdot[H2O]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*q_nocor;
        J[12] += dqdci;               /* dwdot[H]/d[H2] */
        J[15] += dqdci;               /* dwdot[OH]/d[H2] */
        J[16] -= dqdci;               /* dwdot[H2O]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*sc[0];
        J[36] += dqdci;               /* dwdot[H]/d[OH] */
        J[39] += dqdci;               /* dwdot[OH]/d[OH] */
        J[40] -= dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*q_nocor + k_f;
        J[48] += dqdci;               /* dwdot[H]/d[H2O] */
        J[51] += dqdci;               /* dwdot[OH]/d[H2O] */
        J[52] -= dqdci;               /* dwdot[H2O]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (TB[5][4] - 1)*q_nocor;
        J[60] += dqdci;               /* dwdot[H]/d[O2] */
        J[63] += dqdci;               /* dwdot[OH]/d[O2] */
        J[64] -= dqdci;               /* dwdot[H2O]/d[O2] */
        /* d()/d[N2] */
        dqdci = (TB[5][3] - 1)*q_nocor;
        J[96] += dqdci;               /* dwdot[H]/d[N2] */
        J[99] += dqdci;               /* dwdot[OH]/d[N2] */
        J[100] -= dqdci;              /* dwdot[H2O]/d[N2] */
        /* d()/d[HE] */
        dqdci = (TB[5][2] - 1)*q_nocor;
        J[120] += dqdci;              /* dwdot[H]/d[HE] */
        J[123] += dqdci;              /* dwdot[OH]/d[HE] */
        J[124] -= dqdci;              /* dwdot[H2O]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac - k_r*sc[3];
        dqdc[1] = TB[5][0];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac - k_r*sc[0];
        dqdc[4] = TB[5][1] + k_f;
        dqdc[5] = TB[5][4];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[5][3];
        dqdc[9] = dcdc_fac;
        dqdc[10] = TB[5][2];
        for (int k=0; k<11; k++) {
            J[12*k+0] += dqdc[k];
            J[12*k+3] += dqdc[k];
            J[12*k+4] -= dqdc[k];
        }
    }
    J[132] += dqdT; /* dwdot[H]/dT */
    J[135] += dqdT; /* dwdot[OH]/dT */
    J[136] -= dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: O + OH + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[1] + (TB[6][1] - 1)*sc[4] + (TB[6][2] - 1)*sc[9] + (TB[6][3] - 1)*sc[10];
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[2] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[6][0] - 1)*q_nocor;
        J[14] -= dqdci;               /* dwdot[O]/d[H2] */
        J[15] -= dqdci;               /* dwdot[OH]/d[H2] */
        J[18] += dqdci;               /* dwdot[HO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[3];
        J[26] -= dqdci;               /* dwdot[O]/d[O] */
        J[27] -= dqdci;               /* dwdot[OH]/d[O] */
        J[30] += dqdci;               /* dwdot[HO2]/d[O] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[2];
        J[38] -= dqdci;               /* dwdot[O]/d[OH] */
        J[39] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[42] += dqdci;               /* dwdot[HO2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[6][1] - 1)*q_nocor;
        J[50] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[51] -= dqdci;               /* dwdot[OH]/d[H2O] */
        J[54] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[74] -= dqdci;               /* dwdot[O]/d[HO2] */
        J[75] -= dqdci;               /* dwdot[OH]/d[HO2] */
        J[78] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[AR] */
        dqdci = (TB[6][2] - 1)*q_nocor;
        J[110] -= dqdci;              /* dwdot[O]/d[AR] */
        J[111] -= dqdci;              /* dwdot[OH]/d[AR] */
        J[114] += dqdci;              /* dwdot[HO2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[6][3] - 1)*q_nocor;
        J[122] -= dqdci;              /* dwdot[O]/d[HE] */
        J[123] -= dqdci;              /* dwdot[OH]/d[HE] */
        J[126] += dqdci;              /* dwdot[HO2]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[6][0];
        dqdc[2] = dcdc_fac + k_f*sc[3];
        dqdc[3] = dcdc_fac + k_f*sc[2];
        dqdc[4] = TB[6][1];
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[6][2];
        dqdc[10] = TB[6][3];
        for (int k=0; k<11; k++) {
            J[12*k+2] -= dqdc[k];
            J[12*k+3] -= dqdc[k];
            J[12*k+6] += dqdc[k];
        }
    }
    J[134] -= dqdT; /* dwdot[O]/dT */
    J[135] -= dqdT; /* dwdot[OH]/dT */
    J[138] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 8: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[5] -= dqdci;                /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[3];
    J[24] -= dqdci;               /* dwdot[H]/d[O] */
    J[26] += dqdci;               /* dwdot[O]/d[O] */
    J[27] += dqdci;               /* dwdot[OH]/d[O] */
    J[29] -= dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[36] -= dqdci;               /* dwdot[H]/d[OH] */
    J[38] += dqdci;               /* dwdot[O]/d[OH] */
    J[39] += dqdci;               /* dwdot[OH]/d[OH] */
    J[41] -= dqdci;               /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[60] -= dqdci;               /* dwdot[H]/d[O2] */
    J[62] += dqdci;               /* dwdot[O]/d[O2] */
    J[63] += dqdci;               /* dwdot[OH]/d[O2] */
    J[65] -= dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[132] -= dqdT;               /* dwdot[H]/dT */
    J[134] += dqdT;               /* dwdot[O]/dT */
    J[135] += dqdT;               /* dwdot[OH]/dT */
    J[137] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 9: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[12] += dqdci;               /* dwdot[H]/d[H2] */
    J[13] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[14] -= dqdci;               /* dwdot[O]/d[H2] */
    J[15] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[24] += dqdci;               /* dwdot[H]/d[O] */
    J[25] -= dqdci;               /* dwdot[H2]/d[O] */
    J[26] -= dqdci;               /* dwdot[O]/d[O] */
    J[27] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[36] += dqdci;               /* dwdot[H]/d[OH] */
    J[37] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[38] -= dqdci;               /* dwdot[O]/d[OH] */
    J[39] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[132] += dqdT;               /* dwdot[H]/dT */
    J[133] -= dqdT;               /* dwdot[H2]/dT */
    J[134] -= dqdT;               /* dwdot[O]/dT */
    J[135] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 10: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[12] += dqdci;               /* dwdot[H]/d[H2] */
    J[13] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[14] -= dqdci;               /* dwdot[O]/d[H2] */
    J[15] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[24] += dqdci;               /* dwdot[H]/d[O] */
    J[25] -= dqdci;               /* dwdot[H2]/d[O] */
    J[26] -= dqdci;               /* dwdot[O]/d[O] */
    J[27] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[36] += dqdci;               /* dwdot[H]/d[OH] */
    J[37] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[38] -= dqdci;               /* dwdot[O]/d[OH] */
    J[39] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[132] += dqdT;               /* dwdot[H]/dT */
    J[133] -= dqdT;               /* dwdot[H2]/dT */
    J[134] -= dqdT;               /* dwdot[O]/dT */
    J[135] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[4];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[3] -= dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[12] += dqdci;               /* dwdot[H]/d[H2] */
    J[13] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[15] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[16] += dqdci;               /* dwdot[H2O]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[36] += dqdci;               /* dwdot[H]/d[OH] */
    J[37] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[39] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[40] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0];
    J[48] += dqdci;               /* dwdot[H]/d[H2O] */
    J[49] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[51] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[52] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[132] += dqdT;               /* dwdot[H]/dT */
    J[133] -= dqdT;               /* dwdot[H2]/dT */
    J[135] -= dqdT;               /* dwdot[OH]/dT */
    J[136] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 12: OH + OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[3];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= 2 * q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[26] += dqdci;               /* dwdot[O]/d[O] */
    J[27] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[28] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2*sc[3];
    J[38] += dqdci;               /* dwdot[O]/d[OH] */
    J[39] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[40] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[50] += dqdci;               /* dwdot[O]/d[H2O] */
    J[51] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[52] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[134] += dqdT;               /* dwdot[O]/dT */
    J[135] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[136] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 13: H2 + AR <=> H + H + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[0]*sc[9];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (2*h_RT[0] + h_RT[9]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2*sc[0]*sc[9];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[12] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[13] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1] - k_r*sc[0]*sc[0];
    J[108] += 2 * dqdci;          /* dwdot[H]/d[AR] */
    J[109] -= dqdci;              /* dwdot[H2]/d[AR] */
    /* d()/dT */
    J[132] += 2 * dqdT;           /* dwdot[H]/dT */
    J[133] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 14: H2 + HE <=> H + H + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[0]*sc[10];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (2*h_RT[0] + h_RT[10]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2*sc[0]*sc[10];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[10];
    J[12] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[13] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[HE] */
    dqdci =  + k_f*sc[1] - k_r*sc[0]*sc[0];
    J[120] += 2 * dqdci;          /* dwdot[H]/d[HE] */
    J[121] -= dqdci;              /* dwdot[H2]/d[HE] */
    /* d()/dT */
    J[132] += 2 * dqdT;           /* dwdot[H]/dT */
    J[133] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 15: O + O + AR <=> O2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[2]*sc[9];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[2] + h_RT[9]) + (h_RT[5] + h_RT[9]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2*sc[2]*sc[9];
    J[26] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[29] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[62] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[65] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[2]*sc[2] - k_r*sc[5];
    J[110] += -2 * dqdci;         /* dwdot[O]/d[AR] */
    J[113] += dqdci;              /* dwdot[O2]/d[AR] */
    /* d()/dT */
    J[134] += -2 * dqdT;          /* dwdot[O]/dT */
    J[137] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 16: O + O + HE <=> O2 + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[2]*sc[10];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[10];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[10]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2*sc[2]*sc[10];
    J[26] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[29] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[10];
    J[62] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[65] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[HE] */
    dqdci =  + k_f*sc[2]*sc[2] - k_r*sc[5];
    J[122] += -2 * dqdci;         /* dwdot[O]/d[HE] */
    J[125] += dqdci;              /* dwdot[O2]/d[HE] */
    /* d()/dT */
    J[134] += -2 * dqdT;          /* dwdot[O]/dT */
    J[137] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 17: H2O + H2O <=> H + OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[4];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3]*sc[4];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[4]) + (h_RT[0] + h_RT[3] + h_RT[4]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3]*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] -= dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0]*sc[4];
    J[36] += dqdci;               /* dwdot[H]/d[OH] */
    J[39] += dqdci;               /* dwdot[OH]/d[OH] */
    J[40] -= dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  + k_f*2*sc[4] - k_r*sc[0]*sc[3];
    J[48] += dqdci;               /* dwdot[H]/d[H2O] */
    J[51] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[52] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[132] += dqdT;               /* dwdot[H]/dT */
    J[135] += dqdT;               /* dwdot[OH]/dT */
    J[136] -= dqdT;               /* dwdot[H2O]/dT */

    /*reaction 18: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[5] += dqdci;                /* dwdot[O2]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[12] -= dqdci;               /* dwdot[H]/d[H2] */
    J[13] += dqdci;               /* dwdot[H2]/d[H2] */
    J[17] += dqdci;               /* dwdot[O2]/d[H2] */
    J[18] -= dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[1];
    J[60] -= dqdci;               /* dwdot[H]/d[O2] */
    J[61] += dqdci;               /* dwdot[H2]/d[O2] */
    J[65] += dqdci;               /* dwdot[O2]/d[O2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[72] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[73] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[77] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[132] -= dqdT;               /* dwdot[H]/dT */
    J[133] += dqdT;               /* dwdot[H2]/dT */
    J[137] += dqdT;               /* dwdot[O2]/dT */
    J[138] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 19: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[3];
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (2*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += 2 * dqdci;            /* dwdot[OH]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[3];
    J[36] -= dqdci;               /* dwdot[H]/d[OH] */
    J[39] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[42] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[72] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[75] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[132] -= dqdT;               /* dwdot[H]/dT */
    J[135] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[138] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 20: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[26] -= dqdci;               /* dwdot[O]/d[O] */
    J[27] += dqdci;               /* dwdot[OH]/d[O] */
    J[29] += dqdci;               /* dwdot[O2]/d[O] */
    J[30] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[38] -= dqdci;               /* dwdot[O]/d[OH] */
    J[39] += dqdci;               /* dwdot[OH]/d[OH] */
    J[41] += dqdci;               /* dwdot[O2]/d[OH] */
    J[42] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[3];
    J[62] -= dqdci;               /* dwdot[O]/d[O2] */
    J[63] += dqdci;               /* dwdot[OH]/d[O2] */
    J[65] += dqdci;               /* dwdot[O2]/d[O2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[74] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[75] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[77] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[134] -= dqdT;               /* dwdot[O]/dT */
    J[135] += dqdT;               /* dwdot[OH]/dT */
    J[137] += dqdT;               /* dwdot[O2]/dT */
    J[138] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[39] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[40] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[41] += dqdci;               /* dwdot[O2]/d[OH] */
    J[42] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[51] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[52] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[54] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[63] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[64] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[65] += dqdci;               /* dwdot[O2]/d[O2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[75] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[76] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[77] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[135] -= dqdT;               /* dwdot[OH]/dT */
    J[136] += dqdT;               /* dwdot[H2O]/dT */
    J[137] += dqdT;               /* dwdot[O2]/dT */
    J[138] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[65] += dqdci;               /* dwdot[O2]/d[O2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[77] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[78] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[79] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[89] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[137] += dqdT;               /* dwdot[O2]/dT */
    J[138] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[139] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[65] += dqdci;               /* dwdot[O2]/d[O2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[77] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[78] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[79] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[89] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[137] += dqdT;               /* dwdot[O2]/dT */
    J[138] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[139] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 24: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[36] -= dqdci;               /* dwdot[H]/d[OH] */
    J[39] += dqdci;               /* dwdot[OH]/d[OH] */
    J[40] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[48] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[51] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[52] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[55] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[87] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[132] -= dqdT;               /* dwdot[H]/dT */
    J[135] += dqdT;               /* dwdot[OH]/dT */
    J[136] += dqdT;               /* dwdot[H2O]/dT */
    J[139] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 25: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[6];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[12] -= dqdci;               /* dwdot[H]/d[H2] */
    J[13] += dqdci;               /* dwdot[H2]/d[H2] */
    J[18] += dqdci;               /* dwdot[HO2]/d[H2] */
    J[19] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[1];
    J[72] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[73] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[78] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[79] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[85] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[132] -= dqdT;               /* dwdot[H]/dT */
    J[133] += dqdT;               /* dwdot[H2]/dT */
    J[138] += dqdT;               /* dwdot[HO2]/dT */
    J[139] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[3] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[26] -= dqdci;               /* dwdot[O]/d[O] */
    J[27] += dqdci;               /* dwdot[OH]/d[O] */
    J[30] += dqdci;               /* dwdot[HO2]/d[O] */
    J[31] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[38] -= dqdci;               /* dwdot[O]/d[OH] */
    J[39] += dqdci;               /* dwdot[OH]/d[OH] */
    J[42] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[74] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[75] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[78] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[79] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[86] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[87] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[134] -= dqdT;               /* dwdot[O]/dT */
    J[135] += dqdT;               /* dwdot[OH]/dT */
    J[138] += dqdT;               /* dwdot[HO2]/dT */
    J[139] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[39] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[40] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[42] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[51] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[52] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[54] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[55] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[75] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[76] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[78] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[79] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[87] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[135] -= dqdT;               /* dwdot[OH]/dT */
    J[136] += dqdT;               /* dwdot[H2O]/dT */
    J[138] += dqdT;               /* dwdot[HO2]/dT */
    J[139] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[39] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[40] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[42] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[51] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[52] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[54] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[55] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[75] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[76] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[78] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[79] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[87] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[135] -= dqdT;               /* dwdot[OH]/dT */
    J[136] += dqdT;               /* dwdot[H2O]/dT */
    J[138] += dqdT;               /* dwdot[HO2]/dT */
    J[139] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 29: HO2 + H <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[4] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[O]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[24] -= dqdci;               /* dwdot[H]/d[O] */
    J[26] += dqdci;               /* dwdot[O]/d[O] */
    J[28] += dqdci;               /* dwdot[H2O]/d[O] */
    J[30] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[48] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[50] += dqdci;               /* dwdot[O]/d[H2O] */
    J[52] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[54] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[72] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[74] += dqdci;               /* dwdot[O]/d[HO2] */
    J[76] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[132] -= dqdT;               /* dwdot[H]/dT */
    J[134] += dqdT;               /* dwdot[O]/dT */
    J[136] += dqdT;               /* dwdot[H2O]/dT */
    J[138] -= dqdT;               /* dwdot[HO2]/dT */

    double c_R[11], dcRdT[11], e_RT[11];
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
    for (int k = 0; k < 11; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[132+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 11; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 11; ++m) {
            dehmixdc += eh_RT[m]*J[k*12+m];
        }
        J[k*12+11] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[143] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
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
        /*species 1: H2 */
        species[1] =
            +8.24944200e-04
            -1.62860300e-06 * tc[1]
            -2.84263020e-10 * tc[2]
            +1.65394880e-12 * tc[3];
        /*species 2: O */
        species[2] =
            -1.63816600e-03
            +4.84206400e-06 * tc[1]
            -4.80852900e-09 * tc[2]
            +1.55627840e-12 * tc[3];
        /*species 3: OH */
        species[3] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 4: H2O */
        species[4] =
            +3.47498200e-03
            -1.27093920e-05 * tc[1]
            +2.09057430e-08 * tc[2]
            -1.00263520e-11 * tc[3];
        /*species 5: O2 */
        species[5] =
            +1.12748600e-03
            -1.15123000e-06 * tc[1]
            +3.94163100e-09 * tc[2]
            -3.50742160e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +6.56922600e-03
            -2.97002600e-07 * tc[1]
            -1.38774180e-08 * tc[2]
            +9.88606000e-12 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
        /*species 9: AR */
        species[9] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 10: HE */
        species[10] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
    } else {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: H2 */
        species[1] =
            +7.00064400e-04
            -1.12676580e-07 * tc[1]
            -2.76947340e-11 * tc[2]
            +6.33100800e-15 * tc[3];
        /*species 2: O */
        species[2] =
            -2.75506200e-05
            -6.20560600e-09 * tc[1]
            +1.36532010e-11 * tc[2]
            -1.74722080e-15 * tc[3];
        /*species 3: OH */
        species[3] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 4: H2O */
        species[4] =
            +3.05629300e-03
            -1.74605200e-06 * tc[1]
            +3.60298800e-10 * tc[2]
            -2.55664720e-14 * tc[3];
        /*species 5: O2 */
        species[5] =
            +6.13519700e-04
            -2.51768400e-07 * tc[1]
            +5.32584300e-11 * tc[2]
            -4.54574000e-15 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.33613600e-03
            -2.94937800e-06 * tc[1]
            +7.04671200e-10 * tc[2]
            -5.72661600e-14 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 9: AR */
        species[9] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 10: HE */
        species[10] =
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

    double q_f[29], q_r[29];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 29; ++i) {
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
    kc[0] = 1.0 / (refC) * exp((g_RT[0] + g_RT[5]) - (g_RT[6]));

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    kc[1] = refC * exp((g_RT[7]) - (g_RT[3] + g_RT[3]));

    /*reaction 3: H2 + M <=> H + H + M */
    kc[2] = refC * exp((g_RT[1]) - (g_RT[0] + g_RT[0]));

    /*reaction 4: O + O + M <=> O2 + M */
    kc[3] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2]) - (g_RT[5]));

    /*reaction 5: O + H + M <=> OH + M */
    kc[4] = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[3]));

    /*reaction 6: H2O + M <=> H + OH + M */
    kc[5] = refC * exp((g_RT[4]) - (g_RT[0] + g_RT[3]));

    /*reaction 7: O + OH + M <=> HO2 + M */
    kc[6] = 1.0 / (refC) * exp((g_RT[2] + g_RT[3]) - (g_RT[6]));

    /*reaction 8: H + O2 <=> O + OH */
    kc[7] = exp((g_RT[0] + g_RT[5]) - (g_RT[2] + g_RT[3]));

    /*reaction 9: O + H2 <=> H + OH */
    kc[8] = exp((g_RT[2] + g_RT[1]) - (g_RT[0] + g_RT[3]));

    /*reaction 10: O + H2 <=> H + OH */
    kc[9] = exp((g_RT[2] + g_RT[1]) - (g_RT[0] + g_RT[3]));

    /*reaction 11: H2 + OH <=> H2O + H */
    kc[10] = exp((g_RT[1] + g_RT[3]) - (g_RT[4] + g_RT[0]));

    /*reaction 12: OH + OH <=> O + H2O */
    kc[11] = exp((g_RT[3] + g_RT[3]) - (g_RT[2] + g_RT[4]));

    /*reaction 13: H2 + AR <=> H + H + AR */
    kc[12] = refC * exp((g_RT[1] + g_RT[9]) - (g_RT[0] + g_RT[0] + g_RT[9]));

    /*reaction 14: H2 + HE <=> H + H + HE */
    kc[13] = refC * exp((g_RT[1] + g_RT[10]) - (g_RT[0] + g_RT[0] + g_RT[10]));

    /*reaction 15: O + O + AR <=> O2 + AR */
    kc[14] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2] + g_RT[9]) - (g_RT[5] + g_RT[9]));

    /*reaction 16: O + O + HE <=> O2 + HE */
    kc[15] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2] + g_RT[10]) - (g_RT[5] + g_RT[10]));

    /*reaction 17: H2O + H2O <=> H + OH + H2O */
    kc[16] = refC * exp((g_RT[4] + g_RT[4]) - (g_RT[0] + g_RT[3] + g_RT[4]));

    /*reaction 18: HO2 + H <=> H2 + O2 */
    kc[17] = exp((g_RT[6] + g_RT[0]) - (g_RT[1] + g_RT[5]));

    /*reaction 19: HO2 + H <=> OH + OH */
    kc[18] = exp((g_RT[6] + g_RT[0]) - (g_RT[3] + g_RT[3]));

    /*reaction 20: HO2 + O <=> O2 + OH */
    kc[19] = exp((g_RT[6] + g_RT[2]) - (g_RT[5] + g_RT[3]));

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    kc[20] = exp((g_RT[6] + g_RT[3]) - (g_RT[4] + g_RT[5]));

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    kc[21] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[5]));

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    kc[22] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[5]));

    /*reaction 24: H2O2 + H <=> H2O + OH */
    kc[23] = exp((g_RT[7] + g_RT[0]) - (g_RT[4] + g_RT[3]));

    /*reaction 25: H2O2 + H <=> HO2 + H2 */
    kc[24] = exp((g_RT[7] + g_RT[0]) - (g_RT[6] + g_RT[1]));

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    kc[25] = exp((g_RT[7] + g_RT[2]) - (g_RT[3] + g_RT[6]));

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    kc[26] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[4]));

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    kc[27] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[4]));

    /*reaction 29: HO2 + H <=> O + H2O */
    kc[28] = exp((g_RT[6] + g_RT[0]) - (g_RT[2] + g_RT[4]));

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
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -1.012521000000000e+03 * invT
            +6.592218000000000e+00
            -3.298124000000000e+00 * tc[0]
            -4.124721000000000e-04 * tc[1]
            +1.357169166666667e-07 * tc[2]
            +7.896194999999999e-12 * tc[3]
            -2.067436000000000e-14 * tc[4];
        /*species 2: O */
        species[2] =
            +2.914764000000000e+04 * invT
            -1.756599999999997e-02
            -2.946429000000000e+00 * tc[0]
            +8.190830000000000e-04 * tc[1]
            -4.035053333333333e-07 * tc[2]
            +1.335702500000000e-10 * tc[3]
            -1.945348000000000e-14 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.346309130000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.020811000000000e+04 * invT
            +7.966090000000001e-01
            -3.386842000000000e+00 * tc[0]
            -1.737491000000000e-03 * tc[1]
            +1.059116000000000e-06 * tc[2]
            -5.807150833333333e-10 * tc[3]
            +1.253294000000000e-13 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.005249000000000e+03 * invT
            -2.821802000000000e+00
            -3.212936000000000e+00 * tc[0]
            -5.637430000000000e-04 * tc[1]
            +9.593583333333333e-08 * tc[2]
            -1.094897500000000e-10 * tc[3]
            +4.384277000000000e-14 * tc[4];
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
            -1.766315000000000e+04 * invT
            -3.396609000000000e+00
            -3.388754000000000e+00 * tc[0]
            -3.284613000000000e-03 * tc[1]
            +2.475021666666666e-08 * tc[2]
            +3.854838333333333e-10 * tc[3]
            -1.235757500000000e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.453750000000000e+02 * invT
            -1.866001000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.453750000000000e+02 * invT
            +1.584651200000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.350340000000000e+02 * invT
            +4.346533000000000e+00
            -2.991423000000000e+00 * tc[0]
            -3.500322000000000e-04 * tc[1]
            +9.389715000000000e-09 * tc[2]
            +7.692981666666667e-13 * tc[3]
            -7.913760000000000e-17 * tc[4];
        /*species 2: O */
        species[2] =
            +2.923080000000000e+04 * invT
            -2.378248000000000e+00
            -2.542060000000000e+00 * tc[0]
            +1.377531000000000e-05 * tc[1]
            +5.171338333333333e-10 * tc[2]
            -3.792555833333334e-13 * tc[3]
            +2.184026000000000e-17 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.683628750000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 4: H2O */
        species[4] =
            -2.989921000000000e+04 * invT
            -4.190671000000000e+00
            -2.672146000000000e+00 * tc[0]
            -1.528146500000000e-03 * tc[1]
            +1.455043333333333e-07 * tc[2]
            -1.000830000000000e-11 * tc[3]
            +3.195809000000000e-16 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.233930000000000e+03 * invT
            +5.084119999999999e-01
            -3.697578000000000e+00 * tc[0]
            -3.067598500000000e-04 * tc[1]
            +2.098070000000000e-08 * tc[2]
            -1.479400833333333e-12 * tc[3]
            +5.682175000000001e-17 * tc[4];
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
            -1.800696000000000e+04 * invT
            +4.072030000000000e+00
            -4.573167000000000e+00 * tc[0]
            -2.168068000000000e-03 * tc[1]
            +2.457815000000000e-07 * tc[2]
            -1.957420000000000e-11 * tc[3]
            +7.158270000000000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.453750000000000e+02 * invT
            -1.866001000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.453750000000000e+02 * invT
            +1.584651100000000e+00
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
        /*species 0: H */
        species[0] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -1.01252100e+03 * invT
            +5.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91476400e+04 * invT
            -1.01756600e+00
            -2.94642900e+00 * tc[0]
            +8.19083000e-04 * tc[1]
            -4.03505333e-07 * tc[2]
            +1.33570250e-10 * tc[3]
            -1.94534800e-14 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.34630913e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.02081100e+04 * invT
            -2.03391000e-01
            -3.38684200e+00 * tc[0]
            -1.73749100e-03 * tc[1]
            +1.05911600e-06 * tc[2]
            -5.80715083e-10 * tc[3]
            +1.25329400e-13 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.00524900e+03 * invT
            -3.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
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
            -1.76631500e+04 * invT
            -4.39660900e+00
            -3.38875400e+00 * tc[0]
            -3.28461300e-03 * tc[1]
            +2.47502167e-08 * tc[2]
            +3.85483833e-10 * tc[3]
            -1.23575750e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.45375000e+02 * invT
            -2.86600100e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.45375000e+02 * invT
            +5.84651200e-01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.35034000e+02 * invT
            +3.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92308000e+04 * invT
            -3.37824800e+00
            -2.54206000e+00 * tc[0]
            +1.37753100e-05 * tc[1]
            +5.17133833e-10 * tc[2]
            -3.79255583e-13 * tc[3]
            +2.18402600e-17 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.68362875e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 4: H2O */
        species[4] =
            -2.98992100e+04 * invT
            -5.19067100e+00
            -2.67214600e+00 * tc[0]
            -1.52814650e-03 * tc[1]
            +1.45504333e-07 * tc[2]
            -1.00083000e-11 * tc[3]
            +3.19580900e-16 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.23393000e+03 * invT
            -4.91588000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
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
            -1.80069600e+04 * invT
            +3.07203000e+00
            -4.57316700e+00 * tc[0]
            -2.16806800e-03 * tc[1]
            +2.45781500e-07 * tc[2]
            -1.95742000e-11 * tc[3]
            +7.15827000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.45375000e+02 * invT
            -2.86600100e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.45375000e+02 * invT
            +5.84651100e-01
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
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +2.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 2: O */
        species[2] =
            +1.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +2.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 5: O2 */
        species[5] =
            +2.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
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
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +1.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 2: O */
        species[2] =
            +1.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 3: OH */
        species[3] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 4: H2O */
        species[4] =
            +1.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 5: O2 */
        species[5] =
            +2.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
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
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
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
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 2: O */
        species[2] =
            +2.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 3: OH */
        species[3] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +3.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
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
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 2: O */
        species[2] =
            +2.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 3: OH */
        species[3] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 4: H2O */
        species[4] =
            +2.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
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
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
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
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 2: O */
        species[2] =
            +1.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 3: OH */
        species[3] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +2.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +2.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
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
            +2.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +1.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 2: O */
        species[2] =
            +1.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 3: OH */
        species[3] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +1.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +2.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
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
            +3.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
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
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 2: O */
        species[2] =
            +2.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 3: OH */
        species[3] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +3.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
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
            +3.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 2: O */
        species[2] =
            +2.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 3: OH */
        species[3] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +2.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
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
            +4.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
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
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00 * tc[0]
            +8.24944200e-04 * tc[1]
            -4.07150750e-07 * tc[2]
            -3.15847800e-11 * tc[3]
            +1.03371800e-13 * tc[4]
            -3.29409400e+00 ;
        /*species 2: O */
        species[2] =
            +2.94642900e+00 * tc[0]
            -1.63816600e-03 * tc[1]
            +1.21051600e-06 * tc[2]
            -5.34281000e-10 * tc[3]
            +9.72674000e-14 * tc[4]
            +2.96399500e+00 ;
        /*species 3: OH */
        species[3] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
        /*species 4: H2O */
        species[4] =
            +3.38684200e+00 * tc[0]
            +3.47498200e-03 * tc[1]
            -3.17734800e-06 * tc[2]
            +2.32286033e-09 * tc[3]
            -6.26647000e-13 * tc[4]
            +2.59023300e+00 ;
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00 * tc[0]
            +1.12748600e-03 * tc[1]
            -2.87807500e-07 * tc[2]
            +4.37959000e-10 * tc[3]
            -2.19213850e-13 * tc[4]
            +6.03473800e+00 ;
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
        /*species 9: AR */
        species[9] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600100e+00 ;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +9.15348800e-01 ;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00 * tc[0]
            +7.00064400e-04 * tc[1]
            -2.81691450e-08 * tc[2]
            -3.07719267e-12 * tc[3]
            +3.95688000e-16 * tc[4]
            -1.35511000e+00 ;
        /*species 2: O */
        species[2] =
            +2.54206000e+00 * tc[0]
            -2.75506200e-05 * tc[1]
            -1.55140150e-09 * tc[2]
            +1.51702233e-12 * tc[3]
            -1.09201300e-16 * tc[4]
            +4.92030800e+00 ;
        /*species 3: OH */
        species[3] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
        /*species 4: H2O */
        species[4] =
            +2.67214600e+00 * tc[0]
            +3.05629300e-03 * tc[1]
            -4.36513000e-07 * tc[2]
            +4.00332000e-11 * tc[3]
            -1.59790450e-15 * tc[4]
            +6.86281700e+00 ;
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00 * tc[0]
            +6.13519700e-04 * tc[1]
            -6.29421000e-08 * tc[2]
            +5.91760333e-12 * tc[3]
            -2.84108750e-16 * tc[4]
            +3.18916600e+00 ;
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
        /*species 9: AR */
        species[9] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600100e+00 ;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +9.15348900e-01 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * restrict wt)
{
    wt[0] = 1.007970; /*H */
    wt[1] = 2.015940; /*H2 */
    wt[2] = 15.999400; /*O */
    wt[3] = 17.007370; /*OH */
    wt[4] = 18.015340; /*H2O */
    wt[5] = 31.998800; /*O2 */
    wt[6] = 33.006770; /*HO2 */
    wt[7] = 34.014740; /*H2O2 */
    wt[8] = 28.013400; /*N2 */
    wt[9] = 39.948000; /*AR */
    wt[10] = 4.002600; /*HE */

    return;
}


/*save atomic weights into array */
void atomicWeight(double * restrict awt)
{
    awt[0] = 1.007970; /*H */
    awt[1] = 15.999400; /*O */
    awt[2] = 14.006700; /*N */
    awt[3] = 39.948000; /*AR */
    awt[4] = 4.002600; /*HE */
    awt[5] = 12.011150; /*C */

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
  *LENIMC =           47;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC) {
  *LENRMC =         2728;}
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
  *KK =           11;}
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
  WT[           0] =   0.1007969975471497E+01;
  WT[           1] =   0.2015939950942993E+01;
  WT[           2] =   0.1599940013885498E+02;
  WT[           3] =   0.1700737011432648E+02;
  WT[           4] =   0.1801534008979797E+02;
  WT[           5] =   0.3199880027770996E+02;
  WT[           6] =   0.3300677025318146E+02;
  WT[           7] =   0.3401474022865295E+02;
  WT[           8] =   0.2801339912414551E+02;
  WT[           9] =   0.3994800186157227E+02;
  WT[          10] =   0.4002600193023682E+01;
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
  EPS[           1] =   0.3800000000000000E+02;
  EPS[           2] =   0.8000000000000000E+02;
  EPS[           3] =   0.8000000000000000E+02;
  EPS[           4] =   0.5724000000000000E+03;
  EPS[           5] =   0.1074000000000000E+03;
  EPS[           6] =   0.1074000000000000E+03;
  EPS[           7] =   0.1074000000000000E+03;
  EPS[           8] =   0.9753000000000000E+02;
  EPS[           9] =   0.1365000000000000E+03;
  EPS[          10] =   0.1020000000000000E+02;
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
  SIG[           1] =   0.2920000000000000E+01;
  SIG[           2] =   0.2750000000000000E+01;
  SIG[           3] =   0.2750000000000000E+01;
  SIG[           4] =   0.2605000000000000E+01;
  SIG[           5] =   0.3458000000000000E+01;
  SIG[           6] =   0.3458000000000000E+01;
  SIG[           7] =   0.3458000000000000E+01;
  SIG[           8] =   0.3621000000000000E+01;
  SIG[           9] =   0.3330000000000000E+01;
  SIG[          10] =   0.2576000000000000E+01;
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
  DIP[           4] =   0.1844000000000000E+01;
  DIP[           5] =   0.0000000000000000E+00;
  DIP[           6] =   0.0000000000000000E+00;
  DIP[           7] =   0.0000000000000000E+00;
  DIP[           8] =   0.0000000000000000E+00;
  DIP[           9] =   0.0000000000000000E+00;
  DIP[          10] =   0.0000000000000000E+00;
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
  POL[           1] =   0.7900000000000000E+00;
  POL[           2] =   0.0000000000000000E+00;
  POL[           3] =   0.0000000000000000E+00;
  POL[           4] =   0.0000000000000000E+00;
  POL[           5] =   0.1600000000000000E+01;
  POL[           6] =   0.0000000000000000E+00;
  POL[           7] =   0.0000000000000000E+00;
  POL[           8] =   0.1760000000000000E+01;
  POL[           9] =   0.0000000000000000E+00;
  POL[          10] =   0.0000000000000000E+00;
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
  ZROT[           1] =   0.2800000000000000E+03;
  ZROT[           2] =   0.0000000000000000E+00;
  ZROT[           3] =   0.0000000000000000E+00;
  ZROT[           4] =   0.4000000000000000E+01;
  ZROT[           5] =   0.3800000000000000E+01;
  ZROT[           6] =   0.1000000000000000E+01;
  ZROT[           7] =   0.3800000000000000E+01;
  ZROT[           8] =   0.4000000000000000E+01;
  ZROT[           9] =   0.0000000000000000E+00;
  ZROT[          10] =   0.0000000000000000E+00;
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
  NLIN[           1] =            1;
  NLIN[           2] =            0;
  NLIN[           3] =            1;
  NLIN[           4] =            2;
  NLIN[           5] =            1;
  NLIN[           6] =            2;
  NLIN[           7] =            2;
  NLIN[           8] =            1;
  NLIN[           9] =            0;
  NLIN[          10] =            0;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
  COFLAM[           0] =  -0.3253788579319486E+00;
  COFLAM[           1] =   0.3416397862210875E+01;
  COFLAM[           2] =  -0.3631104489846262E+00;
  COFLAM[           3] =   0.1585986202702918E-01;
  COFLAM[           4] =   0.1109465814083037E+02;
  COFLAM[           5] =  -0.1314981109716795E+01;
  COFLAM[           6] =   0.2434839729070233E+00;
  COFLAM[           7] =  -0.8971808760305891E-02;
  COFLAM[           8] =   0.1969666526185763E+01;
  COFLAM[           9] =   0.1801396425547300E+01;
  COFLAM[          10] =  -0.1549102264837655E+00;
  COFLAM[          11] =   0.6908759721669425E-02;
  COFLAM[          12] =   0.1605449323860262E+02;
  COFLAM[          13] =  -0.4103244048729492E+01;
  COFLAM[          14] =   0.6631513498156665E+00;
  COFLAM[          15] =  -0.2977137989389932E-01;
  COFLAM[          16] =   0.2212504705282149E+02;
  COFLAM[          17] =  -0.8451497583355520E+01;
  COFLAM[          18] =   0.1459336751072498E+01;
  COFLAM[          19] =  -0.7286084756548136E-01;
  COFLAM[          20] =  -0.2513872083482838E+01;
  COFLAM[          21] =   0.3151650942831091E+01;
  COFLAM[          22] =  -0.3099607063507329E+00;
  COFLAM[          23] =   0.1344786044241103E-01;
  COFLAM[          24] =   0.5546401577805573E+00;
  COFLAM[          25] =   0.1591057931808813E+01;
  COFLAM[          26] =  -0.5282455808284543E-01;
  COFLAM[          27] =   0.4072391521895438E-03;
  COFLAM[          28] =   0.1486260321341417E+01;
  COFLAM[          29] =   0.1062274672219855E+01;
  COFLAM[          30] =   0.5716849496805523E-01;
  COFLAM[          31] =  -0.6382557166697954E-02;
  COFLAM[          32] =   0.1154289596065436E+02;
  COFLAM[          33] =  -0.2911524559760463E+01;
  COFLAM[          34] =   0.5546581827235442E+00;
  COFLAM[          35] =  -0.2750103005663747E-01;
  COFLAM[          36] =  -0.2731910414923424E+01;
  COFLAM[          37] =   0.3271507387309670E+01;
  COFLAM[          38] =  -0.3454038971012363E+00;
  COFLAM[          39] =   0.1513951680124476E-01;
  COFLAM[          40] =   0.7017839888449308E+01;
  COFLAM[          41] =   0.2176250912976264E+00;
  COFLAM[          42] =   0.5619507686207265E-01;
  COFLAM[          43] =  -0.2368322714922462E-02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void egtransetCOFETA(double* COFETA) {
  COFETA[           0] =  -0.1987529414716897E+02;
  COFETA[           1] =   0.3416397862210910E+01;
  COFETA[           2] =  -0.3631104489846305E+00;
  COFETA[           3] =   0.1585986202702936E-01;
  COFETA[           4] =  -0.1376710086380533E+02;
  COFETA[           5] =   0.9708665581008417E+00;
  COFETA[           6] =  -0.4534923959308444E-01;
  COFETA[           7] =   0.2096011921433090E-02;
  COFETA[           8] =  -0.1481563591580279E+02;
  COFETA[           9] =   0.1801396425547327E+01;
  COFETA[          10] =  -0.1549102264837694E+00;
  COFETA[          11] =   0.6908759721669611E-02;
  COFETA[          12] =  -0.1478508813778873E+02;
  COFETA[          13] =   0.1801396425547305E+01;
  COFETA[          14] =  -0.1549102264837663E+00;
  COFETA[          15] =   0.6908759721669463E-02;
  COFETA[          16] =  -0.1187800764707307E+02;
  COFETA[          17] =  -0.7882519505482808E+00;
  COFETA[          18] =   0.3341408170058896E+00;
  COFETA[          19] =  -0.1986366361647418E-01;
  COFETA[          20] =  -0.1681080110872507E+02;
  COFETA[          21] =   0.2522528725237932E+01;
  COFETA[          22] =  -0.2490712798240462E+00;
  COFETA[          23] =   0.1100615806447165E-01;
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
  COFETA[          36] =  -0.1860218543214680E+02;
  COFETA[          37] =   0.3271507387309783E+01;
  COFETA[          38] =  -0.3454038971012522E+00;
  COFETA[          39] =   0.1513951680124550E-01;
  COFETA[          40] =  -0.1115306958556883E+02;
  COFETA[          41] =   0.2176250912977139E+00;
  COFETA[          42] =   0.5619507686206007E-01;
  COFETA[          43] =  -0.2368322714921862E-02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void egtransetCOFD(double* COFD) {
  COFD[           0] =  -0.1476537827842222E+02;
  COFD[           1] =   0.4195448036781197E+01;
  COFD[           2] =  -0.3279615676730872E+00;
  COFD[           3] =   0.1412388797372668E-01;
  COFD[           4] =  -0.1168685516750758E+02;
  COFD[           5] =   0.2883726644585007E+01;
  COFD[           6] =  -0.1637777890869061E+00;
  COFD[           7] =   0.7265872474830623E-02;
  COFD[           8] =  -0.1500301598770682E+02;
  COFD[           9] =   0.4131920691035950E+01;
  COFD[          10] =  -0.3275337549968449E+00;
  COFD[          11] =   0.1442760802547180E-01;
  COFD[          12] =  -0.1502025959280452E+02;
  COFD[          13] =   0.4138327950805945E+01;
  COFD[          14] =  -0.3284006010219734E+00;
  COFD[          15] =   0.1446660099762999E-01;
  COFD[          16] =  -0.1695183404426207E+02;
  COFD[          17] =   0.4416203740445888E+01;
  COFD[          18] =  -0.3114890192543383E+00;
  COFD[          19] =   0.1159663710619628E-01;
  COFD[          20] =  -0.1719486110604737E+02;
  COFD[          21] =   0.4862600572193149E+01;
  COFD[          22] =  -0.4215909171603991E+00;
  COFD[          23] =   0.1846373197001472E-01;
  COFD[          24] =  -0.1720445512239349E+02;
  COFD[          25] =   0.4866330968714804E+01;
  COFD[          26] =  -0.4220901924876364E+00;
  COFD[          27] =   0.1848595807162496E-01;
  COFD[          28] =  -0.1721362788263434E+02;
  COFD[          29] =   0.4869900420159440E+01;
  COFD[          30] =  -0.4225679381202458E+00;
  COFD[          31] =   0.1850722613298993E-01;
  COFD[          32] =  -0.1676897875889624E+02;
  COFD[          33] =   0.4677755283143664E+01;
  COFD[          34] =  -0.3974479445519867E+00;
  COFD[          35] =   0.1741111717763365E-01;
  COFD[          36] =  -0.1715574553730659E+02;
  COFD[          37] =   0.4798805397829151E+01;
  COFD[          38] =  -0.4058992577741924E+00;
  COFD[          39] =   0.1747030892545845E-01;
  COFD[          40] =  -0.9735381521569142E+01;
  COFD[          41] =   0.2191202791663548E+01;
  COFD[          42] =  -0.7485252920875622E-01;
  COFD[          43] =   0.3470152548075659E-02;
  COFD[          44] =  -0.1168685516750758E+02;
  COFD[          45] =   0.2883726644585007E+01;
  COFD[          46] =  -0.1637777890869061E+00;
  COFD[          47] =   0.7265872474830623E-02;
  COFD[          48] =  -0.1023073719095340E+02;
  COFD[          49] =   0.2153598714195162E+01;
  COFD[          50] =  -0.6969019063693850E-01;
  COFD[          51] =   0.3233961733522977E-02;
  COFD[          52] =  -0.1060109579664357E+02;
  COFD[          53] =   0.2157129122185176E+01;
  COFD[          54] =  -0.6524729611183182E-01;
  COFD[          55] =   0.2809597810849289E-02;
  COFD[          56] =  -0.1060412419807472E+02;
  COFD[          57] =   0.2157123706881953E+01;
  COFD[          58] =  -0.6524690383610432E-01;
  COFD[          59] =   0.2809593789707914E-02;
  COFD[          60] =  -0.1787213075631945E+02;
  COFD[          61] =   0.4905112324365878E+01;
  COFD[          62] =  -0.4173872143518051E+00;
  COFD[          63] =   0.1788845489573668E-01;
  COFD[          64] =  -0.1228690337278460E+02;
  COFD[          65] =   0.2739817607547901E+01;
  COFD[          66] =  -0.1455907617836413E+00;
  COFD[          67] =   0.6496697467578658E-02;
  COFD[          68] =  -0.1229066279603360E+02;
  COFD[          69] =   0.2741057082889793E+01;
  COFD[          70] =  -0.1457627022193994E+00;
  COFD[          71] =   0.6504615244946633E-02;
  COFD[          72] =  -0.1229422250798061E+02;
  COFD[          73] =   0.2742232459895579E+01;
  COFD[          74] =  -0.1459257489141948E+00;
  COFD[          75] =   0.6512123400236565E-02;
  COFD[          76] =  -0.1217511489930491E+02;
  COFD[          77] =   0.2679286720576265E+01;
  COFD[          78] =  -0.1373997446113497E+00;
  COFD[          79] =   0.6125379451516252E-02;
  COFD[          80] =  -0.1272604230538407E+02;
  COFD[          81] =   0.2921765441423414E+01;
  COFD[          82] =  -0.1699072228926505E+00;
  COFD[          83] =   0.7582709264066988E-02;
  COFD[          84] =  -0.9853589402094190E+01;
  COFD[          85] =   0.2052962092552261E+01;
  COFD[          86] =  -0.5895182333307221E-01;
  COFD[          87] =   0.2891951845562504E-02;
  COFD[          88] =  -0.1500301598770682E+02;
  COFD[          89] =   0.4131920691035950E+01;
  COFD[          90] =  -0.3275337549968449E+00;
  COFD[          91] =   0.1442760802547180E-01;
  COFD[          92] =  -0.1060109579664357E+02;
  COFD[          93] =   0.2157129122185176E+01;
  COFD[          94] =  -0.6524729611183182E-01;
  COFD[          95] =   0.2809597810849289E-02;
  COFD[          96] =  -0.1329503778398344E+02;
  COFD[          97] =   0.2938989816956702E+01;
  COFD[          98] =  -0.1706012465018490E+00;
  COFD[          99] =   0.7547358484820122E-02;
  COFD[         100] =  -0.1331106228772027E+02;
  COFD[         101] =   0.2939408772392007E+01;
  COFD[         102] =  -0.1706589698617190E+00;
  COFD[         103] =   0.7549998035548285E-02;
  COFD[         104] =  -0.1894546586601941E+02;
  COFD[         105] =   0.4935720676293722E+01;
  COFD[         106] =  -0.4114057447467022E+00;
  COFD[         107] =   0.1723315986098057E-01;
  COFD[         108] =  -0.1474363832488901E+02;
  COFD[         109] =   0.3350090901251298E+01;
  COFD[         110] =  -0.2245382243795846E+00;
  COFD[         111] =   0.9906528312288899E-02;
  COFD[         112] =  -0.1476396639560408E+02;
  COFD[         113] =   0.3356474034011327E+01;
  COFD[         114] =  -0.2254105976116209E+00;
  COFD[         115] =   0.9946130954817655E-02;
  COFD[         116] =  -0.1478376315874751E+02;
  COFD[         117] =   0.3362741024796286E+01;
  COFD[         118] =  -0.2262670743508222E+00;
  COFD[         119] =   0.9985011202484037E-02;
  COFD[         120] =  -0.1442662508885140E+02;
  COFD[         121] =   0.3216422191096490E+01;
  COFD[         122] =  -0.2069561315576079E+00;
  COFD[         123] =   0.9135288598631700E-02;
  COFD[         124] =  -0.1555067129792547E+02;
  COFD[         125] =   0.3655895323186071E+01;
  COFD[         126] =  -0.2636354398546018E+00;
  COFD[         127] =   0.1157316255067020E-01;
  COFD[         128] =  -0.9681129070187980E+01;
  COFD[         129] =   0.1775605691412164E+01;
  COFD[         130] =  -0.1628295103092079E-01;
  COFD[         131] =   0.7262086859513526E-03;
  COFD[         132] =  -0.1502025959280452E+02;
  COFD[         133] =   0.4138327950805945E+01;
  COFD[         134] =  -0.3284006010219734E+00;
  COFD[         135] =   0.1446660099762999E-01;
  COFD[         136] =  -0.1060412419807472E+02;
  COFD[         137] =   0.2157123706881953E+01;
  COFD[         138] =  -0.6524690383610432E-01;
  COFD[         139] =   0.2809593789707914E-02;
  COFD[         140] =  -0.1331106228772027E+02;
  COFD[         141] =   0.2939408772392007E+01;
  COFD[         142] =  -0.1706589698617190E+00;
  COFD[         143] =   0.7549998035548285E-02;
  COFD[         144] =  -0.1332558556199739E+02;
  COFD[         145] =   0.2938989816956671E+01;
  COFD[         146] =  -0.1706012465018444E+00;
  COFD[         147] =   0.7547358484819904E-02;
  COFD[         148] =  -0.1896004007077217E+02;
  COFD[         149] =   0.4935377774711531E+01;
  COFD[         150] =  -0.4113926173034602E+00;
  COFD[         151] =   0.1723402811988101E-01;
  COFD[         152] =  -0.1473449002110577E+02;
  COFD[         153] =   0.3337793989426398E+01;
  COFD[         154] =  -0.2228575541864579E+00;
  COFD[         155] =   0.9830229599517820E-02;
  COFD[         156] =  -0.1475457475365462E+02;
  COFD[         157] =   0.3343986593397501E+01;
  COFD[         158] =  -0.2237039347662408E+00;
  COFD[         159] =   0.9868653787965939E-02;
  COFD[         160] =  -0.1477418610290300E+02;
  COFD[         161] =   0.3350090901251285E+01;
  COFD[         162] =  -0.2245382243795827E+00;
  COFD[         163] =   0.9906528312288795E-02;
  COFD[         164] =  -0.1442063092189782E+02;
  COFD[         165] =   0.3205844724257626E+01;
  COFD[         166] =  -0.2055144159347264E+00;
  COFD[         167] =   0.9070008933822527E-02;
  COFD[         168] =  -0.1553803545411258E+02;
  COFD[         169] =   0.3641692220020330E+01;
  COFD[         170] =  -0.2617120813083521E+00;
  COFD[         171] =   0.1148659824340386E-01;
  COFD[         172] =  -0.9686347735719096E+01;
  COFD[         173] =   0.1775508180457379E+01;
  COFD[         174] =  -0.1627025530906716E-01;
  COFD[         175] =   0.7256647582828849E-03;
  COFD[         176] =  -0.1695183404426207E+02;
  COFD[         177] =   0.4416203740445888E+01;
  COFD[         178] =  -0.3114890192543383E+00;
  COFD[         179] =   0.1159663710619628E-01;
  COFD[         180] =  -0.1787213075631945E+02;
  COFD[         181] =   0.4905112324365878E+01;
  COFD[         182] =  -0.4173872143518051E+00;
  COFD[         183] =   0.1788845489573668E-01;
  COFD[         184] =  -0.1894546586601941E+02;
  COFD[         185] =   0.4935720676293722E+01;
  COFD[         186] =  -0.4114057447467022E+00;
  COFD[         187] =   0.1723315986098057E-01;
  COFD[         188] =  -0.1896004007077217E+02;
  COFD[         189] =   0.4935377774711531E+01;
  COFD[         190] =  -0.4113926173034602E+00;
  COFD[         191] =   0.1723402811988101E-01;
  COFD[         192] =  -0.1301206458003408E+02;
  COFD[         193] =   0.1429168452568504E+01;
  COFD[         194] =   0.1661557715508861E+00;
  COFD[         195] =  -0.1214321823404827E-01;
  COFD[         196] =  -0.2036133619470493E+02;
  COFD[         197] =   0.5195864695910879E+01;
  COFD[         198] =  -0.4301216528920454E+00;
  COFD[         199] =   0.1744936825492251E-01;
  COFD[         200] =  -0.1973436691750415E+02;
  COFD[         201] =   0.4993125184848788E+01;
  COFD[         202] =  -0.4088531920998837E+00;
  COFD[         203] =   0.1672325900844039E-01;
  COFD[         204] =  -0.1972379377029397E+02;
  COFD[         205] =   0.4985700637038792E+01;
  COFD[         206] =  -0.4077156220815392E+00;
  COFD[         207] =   0.1666668649763390E-01;
  COFD[         208] =  -0.2025975101746627E+02;
  COFD[         209] =   0.5176212790334279E+01;
  COFD[         210] =  -0.4308686680573706E+00;
  COFD[         211] =   0.1761066268522406E-01;
  COFD[         212] =  -0.2001704718429840E+02;
  COFD[         213] =   0.5030685873326798E+01;
  COFD[         214] =  -0.4055157556084015E+00;
  COFD[         215] =   0.1625731432523394E-01;
  COFD[         216] =  -0.1286980431258231E+02;
  COFD[         217] =   0.3009638621722941E+01;
  COFD[         218] =  -0.1808900918568948E+00;
  COFD[         219] =   0.8039983754952261E-02;
  COFD[         220] =  -0.1719486110604737E+02;
  COFD[         221] =   0.4862600572193149E+01;
  COFD[         222] =  -0.4215909171603991E+00;
  COFD[         223] =   0.1846373197001472E-01;
  COFD[         224] =  -0.1228690337278460E+02;
  COFD[         225] =   0.2739817607547901E+01;
  COFD[         226] =  -0.1455907617836413E+00;
  COFD[         227] =   0.6496697467578658E-02;
  COFD[         228] =  -0.1474363832488901E+02;
  COFD[         229] =   0.3350090901251298E+01;
  COFD[         230] =  -0.2245382243795846E+00;
  COFD[         231] =   0.9906528312288899E-02;
  COFD[         232] =  -0.1473449002110577E+02;
  COFD[         233] =   0.3337793989426398E+01;
  COFD[         234] =  -0.2228575541864579E+00;
  COFD[         235] =   0.9830229599517820E-02;
  COFD[         236] =  -0.2036133619470493E+02;
  COFD[         237] =   0.5195864695910879E+01;
  COFD[         238] =  -0.4301216528920454E+00;
  COFD[         239] =   0.1744936825492251E-01;
  COFD[         240] =  -0.1579169675646239E+02;
  COFD[         241] =   0.3572143437285479E+01;
  COFD[         242] =  -0.2518469828462104E+00;
  COFD[         243] =   0.1102533331592793E-01;
  COFD[         244] =  -0.1579979030842146E+02;
  COFD[         245] =   0.3572309323030401E+01;
  COFD[         246] =  -0.2518694694768392E+00;
  COFD[         247] =   0.1102634618303224E-01;
  COFD[         248] =  -0.1580828869487550E+02;
  COFD[         249] =   0.3572786632028933E+01;
  COFD[         250] =  -0.2519341709914386E+00;
  COFD[         251] =   0.1102926053643803E-01;
  COFD[         252] =  -0.1544409203507008E+02;
  COFD[         253] =   0.3434913447661099E+01;
  COFD[         254] =  -0.2339977102148624E+00;
  COFD[         255] =   0.1025033359000777E-01;
  COFD[         256] =  -0.1647508018746298E+02;
  COFD[         257] =   0.3812741716466224E+01;
  COFD[         258] =  -0.2815324383266344E+00;
  COFD[         259] =   0.1224615497227191E-01;
  COFD[         260] =  -0.1028043354360960E+02;
  COFD[         261] =   0.1896303471491212E+01;
  COFD[         262] =  -0.3330438429816188E-01;
  COFD[         263] =   0.1524890767209287E-02;
  COFD[         264] =  -0.1720445512239349E+02;
  COFD[         265] =   0.4866330968714804E+01;
  COFD[         266] =  -0.4220901924876364E+00;
  COFD[         267] =   0.1848595807162496E-01;
  COFD[         268] =  -0.1229066279603360E+02;
  COFD[         269] =   0.2741057082889793E+01;
  COFD[         270] =  -0.1457627022193994E+00;
  COFD[         271] =   0.6504615244946633E-02;
  COFD[         272] =  -0.1476396639560408E+02;
  COFD[         273] =   0.3356474034011327E+01;
  COFD[         274] =  -0.2254105976116209E+00;
  COFD[         275] =   0.9946130954817655E-02;
  COFD[         276] =  -0.1475457475365462E+02;
  COFD[         277] =   0.3343986593397501E+01;
  COFD[         278] =  -0.2237039347662408E+00;
  COFD[         279] =   0.9868653787965939E-02;
  COFD[         280] =  -0.1973436691750415E+02;
  COFD[         281] =   0.4993125184848788E+01;
  COFD[         282] =  -0.4088531920998837E+00;
  COFD[         283] =   0.1672325900844039E-01;
  COFD[         284] =  -0.1579979030842146E+02;
  COFD[         285] =   0.3572309323030401E+01;
  COFD[         286] =  -0.2518694694768392E+00;
  COFD[         287] =   0.1102634618303224E-01;
  COFD[         288] =  -0.1580720390088048E+02;
  COFD[         289] =   0.3572143437285478E+01;
  COFD[         290] =  -0.2518469828462103E+00;
  COFD[         291] =   0.1102533331592793E-01;
  COFD[         292] =  -0.1581504405580379E+02;
  COFD[         293] =   0.3572299494956542E+01;
  COFD[         294] =  -0.2518681372333372E+00;
  COFD[         295] =   0.1102628617469197E-01;
  COFD[         296] =  -0.1545394137084459E+02;
  COFD[         297] =   0.3436021618273004E+01;
  COFD[         298] =  -0.2341477712817018E+00;
  COFD[         299] =   0.1025708506984392E-01;
  COFD[         300] =  -0.1647873164232925E+02;
  COFD[         301] =   0.3810732050841141E+01;
  COFD[         302] =  -0.2812626132670425E+00;
  COFD[         303] =   0.1223411282905946E-01;
  COFD[         304] =  -0.1028179959593585E+02;
  COFD[         305] =   0.1896257251270889E+01;
  COFD[         306] =  -0.3329823450147470E-01;
  COFD[         307] =   0.1524623974194515E-02;
  COFD[         308] =  -0.1721362788263434E+02;
  COFD[         309] =   0.4869900420159440E+01;
  COFD[         310] =  -0.4225679381202458E+00;
  COFD[         311] =   0.1850722613298993E-01;
  COFD[         312] =  -0.1229422250798061E+02;
  COFD[         313] =   0.2742232459895579E+01;
  COFD[         314] =  -0.1459257489141948E+00;
  COFD[         315] =   0.6512123400236565E-02;
  COFD[         316] =  -0.1478376315874751E+02;
  COFD[         317] =   0.3362741024796286E+01;
  COFD[         318] =  -0.2262670743508222E+00;
  COFD[         319] =   0.9985011202484037E-02;
  COFD[         320] =  -0.1477418610290300E+02;
  COFD[         321] =   0.3350090901251285E+01;
  COFD[         322] =  -0.2245382243795827E+00;
  COFD[         323] =   0.9906528312288795E-02;
  COFD[         324] =  -0.1972379377029397E+02;
  COFD[         325] =   0.4985700637038792E+01;
  COFD[         326] =  -0.4077156220815392E+00;
  COFD[         327] =   0.1666668649763390E-01;
  COFD[         328] =  -0.1580828869487550E+02;
  COFD[         329] =   0.3572786632028933E+01;
  COFD[         330] =  -0.2519341709914386E+00;
  COFD[         331] =   0.1102926053643803E-01;
  COFD[         332] =  -0.1581504405580379E+02;
  COFD[         333] =   0.3572299494956542E+01;
  COFD[         334] =  -0.2518681372333372E+00;
  COFD[         335] =   0.1102628617469197E-01;
  COFD[         336] =  -0.1582224453447645E+02;
  COFD[         337] =   0.3572143437285498E+01;
  COFD[         338] =  -0.2518469828462133E+00;
  COFD[         339] =   0.1102533331592807E-01;
  COFD[         340] =  -0.1546403419937241E+02;
  COFD[         341] =   0.3437367417375052E+01;
  COFD[         342] =  -0.2343299630280416E+00;
  COFD[         343] =   0.1026528034463288E-01;
  COFD[         344] =  -0.1648290476852054E+02;
  COFD[         345] =   0.3809087920338710E+01;
  COFD[         346] =  -0.2810417778226483E+00;
  COFD[         347] =   0.1222425349132379E-01;
  COFD[         348] =  -0.1028308834746900E+02;
  COFD[         349] =   0.1896214006095980E+01;
  COFD[         350] =  -0.3329247904451228E-01;
  COFD[         351] =   0.1524374377170872E-02;
  COFD[         352] =  -0.1676897875889624E+02;
  COFD[         353] =   0.4677755283143664E+01;
  COFD[         354] =  -0.3974479445519867E+00;
  COFD[         355] =   0.1741111717763365E-01;
  COFD[         356] =  -0.1217511489930491E+02;
  COFD[         357] =   0.2679286720576265E+01;
  COFD[         358] =  -0.1373997446113497E+00;
  COFD[         359] =   0.6125379451516252E-02;
  COFD[         360] =  -0.1442662508885140E+02;
  COFD[         361] =   0.3216422191096490E+01;
  COFD[         362] =  -0.2069561315576079E+00;
  COFD[         363] =   0.9135288598631700E-02;
  COFD[         364] =  -0.1442063092189782E+02;
  COFD[         365] =   0.3205844724257626E+01;
  COFD[         366] =  -0.2055144159347264E+00;
  COFD[         367] =   0.9070008933822527E-02;
  COFD[         368] =  -0.2025975101746627E+02;
  COFD[         369] =   0.5176212790334279E+01;
  COFD[         370] =  -0.4308686680573706E+00;
  COFD[         371] =   0.1761066268522406E-01;
  COFD[         372] =  -0.1544409203507008E+02;
  COFD[         373] =   0.3434913447661099E+01;
  COFD[         374] =  -0.2339977102148624E+00;
  COFD[         375] =   0.1025033359000777E-01;
  COFD[         376] =  -0.1545394137084459E+02;
  COFD[         377] =   0.3436021618273004E+01;
  COFD[         378] =  -0.2341477712817018E+00;
  COFD[         379] =   0.1025708506984392E-01;
  COFD[         380] =  -0.1546403419937241E+02;
  COFD[         381] =   0.3437367417375052E+01;
  COFD[         382] =  -0.2343299630280416E+00;
  COFD[         383] =   0.1026528034463288E-01;
  COFD[         384] =  -0.1521415383275665E+02;
  COFD[         385] =   0.3348053449783496E+01;
  COFD[         386] =  -0.2233657260417595E+00;
  COFD[         387] =   0.9817787279109837E-02;
  COFD[         388] =  -0.1621324217553288E+02;
  COFD[         389] =   0.3714519551995751E+01;
  COFD[         390] =  -0.2693102767279618E+00;
  COFD[         391] =   0.1173858202719274E-01;
  COFD[         392] =  -0.1014817848789097E+02;
  COFD[         393] =   0.1823570562174538E+01;
  COFD[         394] =  -0.2296236173006269E-01;
  COFD[         395] =   0.1036389912916652E-02;
  COFD[         396] =  -0.1715574553730659E+02;
  COFD[         397] =   0.4798805397829151E+01;
  COFD[         398] =  -0.4058992577741924E+00;
  COFD[         399] =   0.1747030892545845E-01;
  COFD[         400] =  -0.1272604230538407E+02;
  COFD[         401] =   0.2921765441423414E+01;
  COFD[         402] =  -0.1699072228926505E+00;
  COFD[         403] =   0.7582709264066988E-02;
  COFD[         404] =  -0.1555067129792547E+02;
  COFD[         405] =   0.3655895323186071E+01;
  COFD[         406] =  -0.2636354398546018E+00;
  COFD[         407] =   0.1157316255067020E-01;
  COFD[         408] =  -0.1553803545411258E+02;
  COFD[         409] =   0.3641692220020330E+01;
  COFD[         410] =  -0.2617120813083521E+00;
  COFD[         411] =   0.1148659824340386E-01;
  COFD[         412] =  -0.2001704718429840E+02;
  COFD[         413] =   0.5030685873326798E+01;
  COFD[         414] =  -0.4055157556084015E+00;
  COFD[         415] =   0.1625731432523394E-01;
  COFD[         416] =  -0.1647508018746298E+02;
  COFD[         417] =   0.3812741716466224E+01;
  COFD[         418] =  -0.2815324383266344E+00;
  COFD[         419] =   0.1224615497227191E-01;
  COFD[         420] =  -0.1647873164232925E+02;
  COFD[         421] =   0.3810732050841141E+01;
  COFD[         422] =  -0.2812626132670425E+00;
  COFD[         423] =   0.1223411282905946E-01;
  COFD[         424] =  -0.1648290476852054E+02;
  COFD[         425] =   0.3809087920338710E+01;
  COFD[         426] =  -0.2810417778226483E+00;
  COFD[         427] =   0.1222425349132379E-01;
  COFD[         428] =  -0.1621324217553288E+02;
  COFD[         429] =   0.3714519551995751E+01;
  COFD[         430] =  -0.2693102767279618E+00;
  COFD[         431] =   0.1173858202719274E-01;
  COFD[         432] =  -0.1736905261068727E+02;
  COFD[         433] =   0.4134132180231587E+01;
  COFD[         434] =  -0.3219126846923269E+00;
  COFD[         435] =   0.1394110637546756E-01;
  COFD[         436] =  -0.1075685106039131E+02;
  COFD[         437] =   0.2105974066569993E+01;
  COFD[         438] =  -0.6310245347805579E-01;
  COFD[         439] =   0.2931474019944176E-02;
  COFD[         440] =  -0.9735381521569142E+01;
  COFD[         441] =   0.2191202791663548E+01;
  COFD[         442] =  -0.7485252920875622E-01;
  COFD[         443] =   0.3470152548075659E-02;
  COFD[         444] =  -0.9853589402094190E+01;
  COFD[         445] =   0.2052962092552261E+01;
  COFD[         446] =  -0.5895182333307221E-01;
  COFD[         447] =   0.2891951845562504E-02;
  COFD[         448] =  -0.9681129070187980E+01;
  COFD[         449] =   0.1775605691412164E+01;
  COFD[         450] =  -0.1628295103092079E-01;
  COFD[         451] =   0.7262086859513526E-03;
  COFD[         452] =  -0.9686347735719096E+01;
  COFD[         453] =   0.1775508180457379E+01;
  COFD[         454] =  -0.1627025530906716E-01;
  COFD[         455] =   0.7256647582828849E-03;
  COFD[         456] =  -0.1286980431258231E+02;
  COFD[         457] =   0.3009638621722941E+01;
  COFD[         458] =  -0.1808900918568948E+00;
  COFD[         459] =   0.8039983754952261E-02;
  COFD[         460] =  -0.1028043354360960E+02;
  COFD[         461] =   0.1896303471491212E+01;
  COFD[         462] =  -0.3330438429816188E-01;
  COFD[         463] =   0.1524890767209287E-02;
  COFD[         464] =  -0.1028179959593585E+02;
  COFD[         465] =   0.1896257251270889E+01;
  COFD[         466] =  -0.3329823450147470E-01;
  COFD[         467] =   0.1524623974194515E-02;
  COFD[         468] =  -0.1028308834746900E+02;
  COFD[         469] =   0.1896214006095980E+01;
  COFD[         470] =  -0.3329247904451228E-01;
  COFD[         471] =   0.1524374377170872E-02;
  COFD[         472] =  -0.1014817848789097E+02;
  COFD[         473] =   0.1823570562174538E+01;
  COFD[         474] =  -0.2296236173006269E-01;
  COFD[         475] =   0.1036389912916652E-02;
  COFD[         476] =  -0.1075685106039131E+02;
  COFD[         477] =   0.2105974066569993E+01;
  COFD[         478] =  -0.6310245347805579E-01;
  COFD[         479] =   0.2931474019944176E-02;
  COFD[         480] =  -0.7722235809097614E+01;
  COFD[         481] =   0.1138839639005951E+01;
  COFD[         482] =   0.7226014852152088E-01;
  COFD[         483] =  -0.3322612720077926E-02;
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
  KTDIF[           2] =           11;
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
  COFTD[           4] =   0.1525347878079855E+00;
  COFTD[           5] =   0.5464040665780899E-04;
  COFTD[           6] =  -0.2934125089337968E-07;
  COFTD[           7] =   0.4870919731196246E-11;
  COFTD[           8] =   0.2700102625474920E+00;
  COFTD[           9] =   0.3615551214404962E-03;
  COFTD[          10] =  -0.1807447676772846E-06;
  COFTD[          11] =   0.2753212716087208E-10;
  COFTD[          12] =   0.2720417770528547E+00;
  COFTD[          13] =   0.3642754049836650E-03;
  COFTD[          14] =  -0.1821046627192010E-06;
  COFTD[          15] =   0.2773927453061686E-10;
  COFTD[          16] =  -0.1418836567013430E+00;
  COFTD[          17] =   0.7665588601033144E-03;
  COFTD[          18] =  -0.3065500230258140E-06;
  COFTD[          19] =   0.4029595271283201E-10;
  COFTD[          20] =   0.2204829568680314E+00;
  COFTD[          21] =   0.4801643237618346E-03;
  COFTD[          22] =  -0.2329279629943256E-06;
  COFTD[          23] =   0.3464704623498651E-10;
  COFTD[          24] =   0.2209079675460805E+00;
  COFTD[          25] =   0.4810899053474402E-03;
  COFTD[          26] =  -0.2333769631025201E-06;
  COFTD[          27] =   0.3471383309607501E-10;
  COFTD[          28] =   0.2213085142120510E+00;
  COFTD[          29] =   0.4819622095914184E-03;
  COFTD[          30] =  -0.2338001183446032E-06;
  COFTD[          31] =   0.3477677564298334E-10;
  COFTD[          32] =   0.2407445350007590E+00;
  COFTD[          33] =   0.4453434831756282E-03;
  COFTD[          34] =  -0.2181738909148746E-06;
  COFTD[          35] =   0.3269585309435789E-10;
  COFTD[          36] =   0.1654293333730542E+00;
  COFTD[          37] =   0.5612389631132157E-03;
  COFTD[          38] =  -0.2656505647884942E-06;
  COFTD[          39] =   0.3882296208775095E-10;
  COFTD[          40] =   0.3407625358054880E+00;
  COFTD[          41] =  -0.4040578770176091E-04;
  COFTD[          42] =   0.3278795395318293E-07;
  COFTD[          43] =  -0.6270938141265901E-11;
  COFTD[          44] =  -0.1525347878079855E+00;
  COFTD[          45] =  -0.5464040665780899E-04;
  COFTD[          46] =   0.2934125089337968E-07;
  COFTD[          47] =  -0.4870919731196246E-11;
  COFTD[          48] =   0.0000000000000000E+00;
  COFTD[          49] =   0.0000000000000000E+00;
  COFTD[          50] =   0.0000000000000000E+00;
  COFTD[          51] =   0.0000000000000000E+00;
  COFTD[          52] =   0.4155834523029710E+00;
  COFTD[          53] =   0.1097383911871691E-04;
  COFTD[          54] =  -0.3960224834087406E-08;
  COFTD[          55] =   0.1144145307552153E-11;
  COFTD[          56] =   0.4219325599835782E+00;
  COFTD[          57] =   0.1114149277732201E-04;
  COFTD[          58] =  -0.4020727469049531E-08;
  COFTD[          59] =   0.1161625074178193E-11;
  COFTD[          60] =   0.6020321200392946E-01;
  COFTD[          61] =   0.5615615962686471E-03;
  COFTD[          62] =  -0.2553727790421687E-06;
  COFTD[          63] =   0.3633898304582181E-10;
  COFTD[          64] =   0.4427392097103151E+00;
  COFTD[          65] =   0.7117708739452678E-04;
  COFTD[          66] =  -0.3847681435586731E-07;
  COFTD[          67] =   0.6863235673780089E-11;
  COFTD[          68] =   0.4444526949057779E+00;
  COFTD[          69] =   0.7145255629999482E-04;
  COFTD[          70] =  -0.3862572696699677E-07;
  COFTD[          71] =   0.6889797705021213E-11;
  COFTD[          72] =   0.4460703094919803E+00;
  COFTD[          73] =   0.7171261254133862E-04;
  COFTD[          74] =  -0.3876630782084383E-07;
  COFTD[          75] =   0.6914873573367535E-11;
  COFTD[          76] =   0.4452620885099525E+00;
  COFTD[          77] =   0.4946972054970773E-04;
  COFTD[          78] =  -0.2630235132337235E-07;
  COFTD[          79] =   0.4903063330400376E-11;
  COFTD[          80] =   0.4225303552740193E+00;
  COFTD[          81] =   0.1320842793755543E-03;
  COFTD[          82] =  -0.7122224283164592E-07;
  COFTD[          83] =   0.1195161063668316E-10;
  COFTD[          84] =   0.1616137301763599E+00;
  COFTD[          85] =   0.4741552334083435E-04;
  COFTD[          86] =  -0.1671152373050399E-07;
  COFTD[          87] =  -0.1889821705016664E-11;
  COFTD[          88] =  -0.3407625358054880E+00;
  COFTD[          89] =   0.4040578770176091E-04;
  COFTD[          90] =  -0.3278795395318293E-07;
  COFTD[          91] =   0.6270938141265901E-11;
  COFTD[          92] =  -0.1616137301763599E+00;
  COFTD[          93] =  -0.4741552334083435E-04;
  COFTD[          94] =   0.1671152373050399E-07;
  COFTD[          95] =   0.1889821705016664E-11;
  COFTD[          96] =   0.3315880233940311E+00;
  COFTD[          97] =  -0.1963882557865052E-04;
  COFTD[          98] =   0.3023888393001825E-07;
  COFTD[          99] =  -0.8449980090305466E-11;
  COFTD[         100] =   0.3422032136897760E+00;
  COFTD[         101] =  -0.2026752702741960E-04;
  COFTD[         102] =   0.3120692705763941E-07;
  COFTD[         103] =  -0.8720490905912508E-11;
  COFTD[         104] =   0.2849835858178697E+00;
  COFTD[         105] =   0.1154600116315077E-03;
  COFTD[         106] =  -0.6171979338326076E-07;
  COFTD[         107] =   0.1015042212845429E-10;
  COFTD[         108] =   0.4402209436618440E+00;
  COFTD[         109] =  -0.4837175963856001E-04;
  COFTD[         110] =   0.4660889007292850E-07;
  COFTD[         111] =  -0.1027684286354751E-10;
  COFTD[         112] =   0.4436492507591140E+00;
  COFTD[         113] =  -0.4874846422125629E-04;
  COFTD[         114] =   0.4697186596249824E-07;
  COFTD[         115] =  -0.1035687579663234E-10;
  COFTD[         116] =   0.4468957655172876E+00;
  COFTD[         117] =  -0.4910519334513472E-04;
  COFTD[         118] =   0.4731559438265286E-07;
  COFTD[         119] =  -0.1043266483507838E-10;
  COFTD[         120] =   0.4220100375895677E+00;
  COFTD[         121] =  -0.4140425208641108E-04;
  COFTD[         122] =   0.4387516156683957E-07;
  COFTD[         123] =  -0.1028602437920772E-10;
  COFTD[         124] =   0.4663152819246087E+00;
  COFTD[         125] =  -0.5601505800988478E-04;
  COFTD[         126] =   0.4659876633602610E-07;
  COFTD[         127] =  -0.9136462914858121E-11;
  COFTD[         128] =   0.0000000000000000E+00;
  COFTD[         129] =   0.0000000000000000E+00;
  COFTD[         130] =   0.0000000000000000E+00;
  COFTD[         131] =   0.0000000000000000E+00;
};




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
!
!                     ----- H2 Kinetic Mechanism -----
!                     -----   Version 6-10-2011  -----
!
! (c) Burke, Chaos, Ju, Dryer, and Klippenstein; Princeton University, 2011.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  HOW TO USE THIS MECHANISM:
!
! (*) Due to limitations of CHEMKIN-II format (specifically, an inability to
!     implement temperature-dependent collision efficiencies in falloff
!     reactions) and the lack of fundamental understanding of the mixing rules
!     for the falloff reactions with the bath gases that have different
!     broadening factors, the present implementation represents a compromise
!     (approximate) formulation.  As a consequence,
!
!     PRIOR TO ITS USE IN THE CALCULATIONS, THIS FILE HAS TO BE MODIFIED.
!     DEPENDING ON WHAT BATH GAS (DILUTANT) IS MOST ABUNDANT IN YOUR SYSTEM
!     (THE PRESENT CHOICES ARE N2, AR, OR HE),  YOU  SHOULD UNCOMMENT THE
!     CORRESPONDING BLOCK FOR THE REACTION H+O2(+M)=HO2(+M), AND COMMENT THE
!     BLOCK FOR OTHER DILUTANT(S).  AS GIVEN, THE MAIN DILUTANT IS SET TO BE N2.
!
!
!  HOW TO REFERENCE THIS MECHANISM:
!
!     M.P. Burke, M. Chaos, Y. Ju, F.L. Dryer, S.J. Klippenstein
!        "Comprehensive H2/O2 Kinetic Model for High-Pressure Combustion,"
!        Int. J. Chem. Kinet. (2011).
!
!  FUTURE REVISIONS/UPDATES MAY BE FOUND ON THE FUELS AND COMBUSTION RESEARCH LABORATORY
!  WEBSITE: < http://www.princeton.edu/mae/people/faculty/dryer/homepage/combustion_lab/ >
!
!
!  HOW TO CONTACT THE AUTHORS:
!
!     Dr. Michael P. Burke
!     R122 Building 200
!     Chemical Sciences and Engineering Division
!     Argonne National Laboratory
!     Argonne, IL 60439
!     Email: mpburke@anl.gov
!
!     Prof. Frederick L. Dryer
!     D-329D Engineering Quadrangle
!     Mechanical and Aerospace Engineering
!     Princeton University
!     Princeton, NJ 08544
!     Phone: 609-258-5206
!     Lab:   609-258-0316
!     FAX:   609-258-1939
!     Email: fldryer@princeton.edu
!
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
!
ELEMENTS
H O N AR HE C
END

SPECIES
H        H2       O        OH
H2O      O2       HO2      H2O2     
N2       AR       HE       
END

!*********************************************************************************

THERMO ALL
0300.00  1000.00  5000.00
H                 120186H   1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 0.02547163E+06-0.04601176E+01 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.02547163E+06-0.04601176E+01                   4
H2                121286H   2               G  0300.00   5000.00  1000.00      1
 0.02991423E+02 0.07000644E-02-0.05633829E-06-0.09231578E-10 0.01582752E-13    2
-0.08350340E+04-0.01355110E+02 0.03298124E+02 0.08249442E-02-0.08143015E-05    3
-0.09475434E-09 0.04134872E-11-0.01012521E+05-0.03294094E+02                   4
O                 120186O   1               G  0300.00   5000.00  1000.00      1
 0.02542060E+02-0.02755062E-03-0.03102803E-07 0.04551067E-10-0.04368052E-14    2
 0.02923080E+06 0.04920308E+02 0.02946429E+02-0.01638166E-01 0.02421032E-04    3
-0.01602843E-07 0.03890696E-11 0.02914764E+06 0.02963995E+02                   4
OH                S 9/01O   1H   1    0    0G   200.000  6000.000 1000.        1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.68362875E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.34630913E+03-6.90432960E-01 4.51532273E+03    4
H2O                20387H   2O   1          G  0300.00   5000.00  1000.00      1
 0.02672146E+02 0.03056293E-01-0.08730260E-05 0.01200996E-08-0.06391618E-13    2
-0.02989921E+06 0.06862817E+02 0.03386842E+02 0.03474982E-01-0.06354696E-04    3
 0.06968581E-07-0.02506588E-10-0.03020811E+06 0.02590233E+02                   4
O2                121386O   2               G  0300.00   5000.00  1000.00      1
 0.03697578E+02 0.06135197E-02-0.01258842E-05 0.01775281E-09-0.01136435E-13    2
-0.01233930E+05 0.03189166E+02 0.03212936E+02 0.01127486E-01-0.05756150E-05    3
 0.01313877E-07-0.08768554E-11-0.01005249E+05 0.06034738E+02                   4
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 1.00021620E+04    4
H2O2              120186H   2O   2          G  0300.00   5000.00  1000.00      1
 0.04573167E+02 0.04336136E-01-0.01474689E-04 0.02348904E-08-0.01431654E-12    2
-0.01800696E+06 0.05011370E+01 0.03388754E+02 0.06569226E-01-0.01485013E-05    3
-0.04625806E-07 0.02471515E-10-0.01766315E+06 0.06785363E+02                   4
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.01487977E-01-0.05684761E-05 0.01009704E-08-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.01408240E-01-0.03963222E-04    3
 0.05641515E-07-0.02444855E-10-0.01020900E+05 0.03950372E+02                   4
AR                120186AR  1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366001E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366001E+02                   4
HE                120186HE  1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.09153489E+01 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.09153488E+01                   4
END

!*********************************************************************************

REACTIONS

!======================
!H2-O2 Chain Reactions
!======================

! Hong et al., Proc. Comb. Inst. 33:309-316 (2011)
H+O2 = O+OH                                 	1.04E+14   0.00  1.5286E+04

! Baulch et al., J. Phys. Chem. Ref. Data, 21:411 (1992)
O+H2 = H+OH						3.818E+12  0.00  7.948E+03
   DUPLICATE
O+H2 = H+OH						8.792E+14  0.00  1.917E+04
   DUPLICATE

! Michael and Sutherland, J. Phys. Chem. 92:3853 (1988)
H2+OH = H2O+H						0.216E+09  1.51  0.343E+04

! Baulch et al., J. Phys. Chem. Ref. Data, 21:411 (1992)
OH+OH = O+H2O						3.34E+04   2.42  -1.93E+03

!============================
!H2-O2 Dissociation Reactions
!============================

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) 
H2+M = H+H+M						4.577E+19 -1.40  1.0438E+05
   H2/2.5/ H2O/12/
   AR/0.0/ HE/0.0/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) 
H2+AR = H+H+AR                              	5.840E+18 -1.10  1.0438E+05
H2+HE = H+H+HE                              	5.840E+18 -1.10  1.0438E+05

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) 
O+O+M = O2+M                                	6.165E+15 -0.50  0.000E+00
   H2/2.5/ H2O/12/
   AR/0.0/ HE/0.0/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) 
O+O+AR = O2+AR                              	1.886E+13  0.00 -1.788E+03
O+O+HE = O2+HE                              	1.886E+13  0.00 -1.788E+03

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) 
O+H+M = OH+M                                	4.714E+18 -1.00  0.000E+00
   H2/2.5/  H2O/12/
   AR/0.75/ HE/0.75/

! Srinivasan and Michael, Int. J. Chem. Kinetic. 38 (2006)
! Rate constant is for Ar with efficiencies from Michael et al., J. Phys. Chem. A, 106 (2002)
H2O+M = H+OH+M                              	6.064E+27 -3.322 1.2079E+05
   H2/3.0/  H2O/0.0/
   HE/1.10/ N2/2.00/
   O2/1.5/
! Efficiencies for CO and CO2 taken from Li et al., Int. J. Chem. Kinet. 36:566-575 (2004)

! Srinivasan and Michael, Int. J. Chem. Kinetic. 38 (2006)
H2O+H2O = H+OH+H2O                          	1.006E+26 -2.44  1.2018E+05

!=================================
! Formation and consumption of HO2
!=================================

! High-pressure limit from Troe, Proc. Comb. Inst. 28:1463-1469 (2000)
! Low-pressure  limit from Michael et al., J. Phys. Chem. A 106:5297-5313
! Centering factors from Fernandes et al., Phys. Chem. Chem. Phys. 10:4313-4321 (2008)
!=================================================================================
! MAIN BATH GAS IS N2 (comment this reaction otherwise)
!
H+O2(+M) = HO2(+M)                          	4.65084E+12  0.44  0.000E+00
   LOW/6.366E+20 -1.72  5.248E+02/
   TROE/0.5  1E-30  1E+30/
   H2/2.0/ H2O/14/ O2/0.78/ AR/0.67/ HE/0.8/
!=================================================================================
! MAIN BATH GAS IS AR OR HE (comment this reaction otherwise)
!
!H+O2(+M) = HO2(+M)                         	4.65084E+12  0.44  0.000E+00
!   LOW/9.042E+19 -1.50  4.922E+02/
!   TROE/0.5 1E-30  1E+30/
!   H2/3.0/ H2O/21/ O2/1.1/ CO/2.7/ CO2/5.4/ HE/1.2/ N2/1.5/
!=================================================================================

! Michael et al., Proc. Comb. Inst. 28:1471 (2000)
!HO2+H = H2+O2                                 	3.659E+06  2.09 -1.451E+03
!Scaled by 0.75
HO2+H = H2+O2                                 	2.750E+06  2.09 -1.451E+03

! Mueller et al., Int. J. Chem. Kinetic. 31:113 (1999) 
HO2+H = OH+OH                               	7.079E+13  0.00  2.950E+02

! Fernandez-Ramos and Varandas, J. Phys. Chem. A 106:4077-4083 (2002)
!HO2+O = O2+OH                               	4.750E+10  1.00 -7.2393E+02
!Scaled by 0.60
HO2+O = O2+OH                               	2.850E+10  1.00 -7.2393E+02

! Keyser, J. Phys. Chem. 92:1193 (1988)
HO2+OH = H2O+O2                             	2.890E+13  0.00 -4.970E+02

!=====================================
!Formation and Consumption of H2O2
!=====================================

! Hippler et al., J. Chem. Phys. 93:1755 (1990)
HO2+HO2 = H2O2+O2                           	4.200E+14  0.00  1.1982E+04
   DUPLICATE
HO2+HO2 = H2O2+O2                           	1.300E+11  0.00 -1.6293E+03
   DUPLICATE

! Troe, Combust. Flame,  158:594-601 (2011)
! Rate constant is for Ar
H2O2(+M) = OH+OH(+M)            			2.00E+12   0.90  4.8749E+04
   LOW/2.49E+24 -2.30 4.8749E+04/
   TROE/0.43 1E-30 1E+30/
   H2O/7.5/ 
   N2/1.5/  O2/1.2/
   HE/0.65/ H2O2/7.7/
! Efficiencies for H2 and CO taken from Li et al., Int. J. Chem. Kinet. 36:566-575 (2004)
   H2/3.7/ 

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+H = H2O+OH                             	2.410E+13  0.00  3.970E+03
H2O2+H = HO2+H2                             	4.820E+13  0.00  7.950E+03
H2O2+O = OH+HO2                             	9.550E+06  2.00  3.970E+03

! Hong et al., J. Phys. Chem. A  114 (2010) 5718–5727
H2O2+OH = HO2+H2O                           	1.740E+12  0.00  3.180E+02
   DUPLICATE
H2O2+OH = HO2+H2O                           	7.590E+13  0.00  7.270E+03
   DUPLICATE

!  JBB added reactions  X1 and X6 from paper

HO2+H = O+H2O                                   3.970E+12  0.00  6.710E+02

O+OH+M = HO2+M                                  8.000E+15  0.00  0.000E+00
   H2/2.0/  H2O/12./ AR/0.7/  HE/0.7/  



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
CH2OCH2            1   252.000     4.760     0.000     0.000     1.000
CH2OCH             1   252.000     4.760     0.000     0.000     1.000
CH3CH2CHO          1   252.000     4.760     0.000     0.000     1.000
                                                                                
C4H                1   357.000     5.180     0.000     0.000     1.000          
C4H2               1   357.000     5.180     0.000     0.000     1.000          
H2C4O              2   357.000     5.180     0.000     0.000     1.000 ! JAM    
C4H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *      
iC4H3              2   357.000     5.180     0.000     0.000     1.000 ! JAM    
nC4H3              2   357.000     5.180     0.000     0.000     1.000 ! JAM    
C4H4               2   357.000     5.180     0.000     0.000     1.000 ! JAM    
iC4H5              2   357.000     5.180     0.000     0.000     1.000 ! JAM    
nC4H5              2   357.000     5.180     0.000     0.000     1.000 ! JAM    
C4H5-2             2   357.000     5.180     0.000     0.000     1.000 !
C4H6               2   357.000     5.180     0.000     0.000     1.000         
C4H6-2             2   357.000     5.180     0.000     0.000     1.000
C4H612             2   357.000     5.180     0.000     0.000     1.000 
CH3CHOCH2          2   357.000     5.180     0.000     0.000     1.000
                                                                                
C5H2               1   357.000     5.180     0.000     0.000     1.000          
C5H3               1   357.000     5.180     0.000     0.000     1.000          
C5H5               1   357.000     5.180     0.000     0.000     1.000          
C5H6               1   357.000     5.180     0.000     0.000     1.000          
lC5H7              1   357.000     5.180     0.000     0.000     1.000
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
nC6H5              2   412.300     5.349     0.000     0.000     1.000  !(JAM)  
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
C6H6               2   464.8       5.29      0.00     10.32      0.000  !  benze
C6H5               2   464.8       5.29      0.00     10.32      0.000  !  benze
C6H5CH3            2   495.3       5.68      0.43     12.30      1.000  !
C6H5C2H3           2   546.2       6.00      0.13     15.00      1.000  !
C6H5CH2            2   495.3       5.68      0.43     12.30      1.000  !
C6H5C2H            2   535.6       5.72      0.77     12.00      1.000  !
A2                 2   630.4       6.18      0.00     16.50      1.000  !
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

C5H5O(2,4)         2   494         5.2       1.6       0.0       1.0
C5H5O(1,2)         2   494         5.2       1.6       0.0       1.0
C5H5O(1,3)         2   494         5.2       1.6       0.0       1.0
C4H5               2   329         5.1       0.0       0.0       1.0
c-C4H5             2   329         5.1       0.0       0.0       1.0
C6H5CO             2   593         5.5       2.8       0.0       1.0
C6H5CHO            2   593         5.47      2.8       0.0       1.0
C6H5C2H5           2   485         5.425     0.4       0.0       1.0
C6H4O2             2   485         5.425     0.4       0.0       1.0
HOC6H4CH3          2   567         5.60      1.6       0.0       1.0
C6H5CH2OH          2   572         5.82      1.7       0.0       1.0
bi-C6H5CH2         2   620         7.24      0.0       0.0       1.0
C5H5OH             2   464.800     5.290     0.000    10.320     0.000  !  as C5H4OH, ZD99
C5H4OH             2   464.800     5.290     0.000    10.320     0.000  !  benze
o-C6H4             2   464.8       5.29      0.00     10.32      0.000  !  benze
C6H5C6H5           2   676.5       6.31      0.00     20.00      1.000  !  biphe
OC6H4CH3           2   567         5.6       1.6       0.0       1.000
C10H8              2   630.4       6.18      0.00     16.50      1.000  !  napht
halene
C6H4CH3            2   495.3       5.68      0.43     12.30      1.000  !

                                                                                
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
