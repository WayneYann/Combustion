#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
from collections import defaultdict

from weaver.mills.CMill import CMill

from pyre.units.pressure import atm
from pyre.units.SI import meter, second, mole, kelvin
from pyre.units.length import cm
from pyre.units.energy import cal, kcal, J, kJ, erg
from pyre.handbook.constants.fundamental import avogadro
from pyre.handbook.constants.fundamental import gas_constant as R

import sys

smallnum = 1e-100
R = 8.31451e7 * erg/mole/kelvin
#Rc = 1.987215583 * cal/mole/kelvin
Rc = 1.98721558317399617591 * cal/mole/kelvin
Patm = 1013250.0
sym  = ""
fsym = "_"

class speciesDb:
    def __init__(self, id, name, mwt):
        self.id = id
        self.symbol = name
        self.weight = mwt
        return


class CPickler(CMill):


    def __init__(self):
        CMill.__init__(self)
        self.species = []
        self.nSpecies = 0
        self.reactionIndex = []
        return


    def _setSpecies(self, mechanism):
        """ For internal use """
        import pyre
        periodic = pyre.handbook.periodicTable()
        
        nSpecies = len(mechanism.species())
        self.species = [ 0.0 for x in range(nSpecies) ]
        
        for species in mechanism.species():
            weight = 0.0 
            for elem, coef in species.composition:
                aw = mechanism.element(elem).weight
                if not aw:
                    aw = periodic.symbol(elem.capitalize()).atomicWeight
                weight += coef * aw

            tempsp = speciesDb(species.id, species.symbol, weight)
            self.species[species.id] = tempsp

        self.nSpecies = nSpecies
        return


    def _statics(self,mechanism):
        self._write()
        self._write(self.line(' Inverse molecular weights'))
        self._write('static const double imw[%d] = {' % (self.nSpecies))
        self._indent()
        for i in range(0,self.nSpecies):
            species = self.species[i]
            text = '1.0 / %f' % (species.weight)
            if (i<self.nSpecies-1):
               text += ',  '
            else:
               text += '};  '
            self._write(text + self.line('%s' % species.symbol))
        self._outdent()
        self._write()

        self._write()
        nReactions = len(mechanism.reaction())
        self._write()
        self._write('static double fwd_A[%d], fwd_beta[%d], fwd_Ea[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('static double low_A[%d], low_beta[%d], low_Ea[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('static double rev_A[%d], rev_beta[%d], rev_Ea[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('static double troe_a[%d],troe_Ts[%d], troe_Tss[%d], troe_Tsss[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions))
        self._write('static double sri_a[%d], sri_b[%d], sri_c[%d], sri_d[%d], sri_e[%d];'
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('static double activation_units[%d], prefactor_units[%d], phase_units[%d];'
                    % (nReactions,nReactions,nReactions))
        self._write('static int is_PD[%d], troe_len[%d], sri_len[%d], nTB[%d], *TBid[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('static double *TB[%d];' 
                    % (nReactions))

        self._write()
        self._write('static double fwd_A_DEF[%d], fwd_beta_DEF[%d], fwd_Ea_DEF[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('static double low_A_DEF[%d], low_beta_DEF[%d], low_Ea_DEF[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('static double rev_A_DEF[%d], rev_beta_DEF[%d], rev_Ea_DEF[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('static double troe_a_DEF[%d],troe_Ts_DEF[%d], troe_Tss_DEF[%d], troe_Tsss_DEF[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions))
        self._write('static double sri_a_DEF[%d], sri_b_DEF[%d], sri_c_DEF[%d], sri_d_DEF[%d], sri_e_DEF[%d];'
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('static double activation_units_DEF[%d], prefactor_units_DEF[%d], phase_units_DEF[%d];'
                    % (nReactions,nReactions,nReactions))
        self._write('static int is_PD_DEF[%d], troe_len_DEF[%d], sri_len_DEF[%d], nTB_DEF[%d], *TBid_DEF[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('static double *TB_DEF[%d];' 
                    % (nReactions))

        # build reverse reaction map
        rmap = {}
        for i, reaction in zip(range(nReactions), mechanism.reaction()):
            rmap[reaction.orig_id-1] = i
        
        self._write('static int rxn_map[%d] = {%s};' % (nReactions, ",".join(str(rmap[x]) for x in range(len(rmap)))))


        self._write('')
        self._write('void GET_REACTION_MAP(int *rmap)')
        self._write('{')
        self._indent()
        self._write('for (int i=0; i<%d; ++i) {' % (nReactions))
        self._indent()        
        self._write('rmap[i] = rxn_map[i];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write()
        self._write("")
        self._write("#include <ReactionData.H>")
        self._write("double* GetParamPtr(int                reaction_id,")
        self._write("                    REACTION_PARAMETER param_id,")
        self._write("                    int                species_id,")
        self._write("                    int                get_default)")
        self._write("{")
        self._write("  double* ret = 0;")
        self._write("  if (reaction_id<0 || reaction_id>=%d) {" % (nReactions))
        self._write("    printf(\"Bad reaction id = %d\",reaction_id);")
        self._write("    abort();")
        self._write("  };")
        self._write("  int mrid = rxn_map[reaction_id];")
        self._write()
        self._write("  if (param_id == THIRD_BODY) {")
        self._write("    if (species_id<0 || species_id>=%d) {" % (self.nSpecies))
        self._write("      printf(\"GetParamPtr: Bad species id = %d\",species_id);")
        self._write("      abort();")
        self._write("    }")
        self._write("    if (get_default) {")
        self._write("      for (int i=0; i<nTB_DEF[mrid]; ++i) {")
        self._write("        if (species_id == TBid_DEF[mrid][i]) {")
        self._write("          ret = &(TB_DEF[mrid][i]);")
        self._write("        }")
        self._write("      }")
        self._write("    }")
        self._write("    else {")
        self._write("      for (int i=0; i<nTB[mrid]; ++i) {")
        self._write("        if (species_id == TBid[mrid][i]) {")
        self._write("          ret = &(TB[mrid][i]);")
        self._write("        }")
        self._write("      }")
        self._write("    }")
        self._write("    if (ret == 0) {")
        self._write("      printf(\"GetParamPtr: No TB for reaction id = %d\",reaction_id);")
        self._write("      abort();")
        self._write("    }")
        self._write("  }")
        self._write("  else {")
        self._write("    if (     param_id == FWD_A)     {ret = (get_default ? &(fwd_A_DEF[mrid]) : &(fwd_A[mrid]));}")
        self._write("      else if (param_id == FWD_BETA)  {ret = (get_default ? &(fwd_beta_DEF[mrid]) : &(fwd_beta[mrid]));}")
        self._write("      else if (param_id == FWD_EA)    {ret = (get_default ? &(fwd_Ea_DEF[mrid]) : &(fwd_Ea[mrid]));}")
        self._write("      else if (param_id == LOW_A)     {ret = (get_default ? &(low_A_DEF[mrid]) : &(low_A[mrid]));}")
        self._write("      else if (param_id == LOW_BETA)  {ret = (get_default ? &(low_beta_DEF[mrid]) : &(low_beta[mrid]));}")
        self._write("      else if (param_id == LOW_EA)    {ret = (get_default ? &(low_Ea_DEF[mrid]) : &(low_Ea[mrid]));}")
        self._write("      else if (param_id == REV_A)     {ret = (get_default ? &(rev_A_DEF[mrid]) : &(rev_A[mrid]));}")
        self._write("      else if (param_id == REV_BETA)  {ret = (get_default ? &(rev_beta_DEF[mrid]) : &(rev_beta[mrid]));}")
        self._write("      else if (param_id == REV_EA)    {ret = (get_default ? &(rev_Ea_DEF[mrid]) : &(rev_Ea[mrid]));}")
        self._write("      else if (param_id == TROE_A)    {ret = (get_default ? &(troe_a_DEF[mrid]) : &(troe_a[mrid]));}")
        self._write("      else if (param_id == TROE_TS)   {ret = (get_default ? &(troe_Ts_DEF[mrid]) : &(troe_Ts[mrid]));}")
        self._write("      else if (param_id == TROE_TSS)  {ret = (get_default ? &(troe_Tss_DEF[mrid]) : &(troe_Tss[mrid]));}")
        self._write("      else if (param_id == TROE_TSSS) {ret = (get_default ? &(troe_Tsss_DEF[mrid]) : &(troe_Tsss[mrid]));}")
        self._write("      else if (param_id == SRI_A)     {ret = (get_default ? &(sri_a_DEF[mrid]) : &(sri_a[mrid]));}")
        self._write("      else if (param_id == SRI_B)     {ret = (get_default ? &(sri_b_DEF[mrid]) : &(sri_b[mrid]));}")
        self._write("      else if (param_id == SRI_C)     {ret = (get_default ? &(sri_c_DEF[mrid]) : &(sri_c[mrid]));}")
        self._write("      else if (param_id == SRI_D)     {ret = (get_default ? &(sri_d_DEF[mrid]) : &(sri_d[mrid]));}")
        self._write("      else if (param_id == SRI_E)     {ret = (get_default ? &(sri_e_DEF[mrid]) : &(sri_e[mrid]));}")
        self._write("    else {")
        self._write("      printf(\"GetParamPtr: Unknown parameter id\");")
        self._write("      abort();")
        self._write("    }")
        self._write("  }")
        self._write("  return ret;")
        self._write("}")
        self._write()
        self._write("void ResetAllParametersToDefault()")
        self._write("{")
        self._write("    for (int i=0; i<%d; i++) {" % (nReactions))
        self._write("        if (nTB[i] != 0) {")
        self._write("            nTB[i] = 0;")
        self._write("            free(TB[i]);")
        self._write("            free(TBid[i]);")
        self._write("        }")
        self._write("")
        self._write("        fwd_A[i]    = fwd_A_DEF[i];")
        self._write("        fwd_beta[i] = fwd_beta_DEF[i];")
        self._write("        fwd_Ea[i]   = fwd_Ea_DEF[i];")
        self._write("")
        self._write("        low_A[i]    = low_A_DEF[i];")
        self._write("        low_beta[i] = low_beta_DEF[i];")
        self._write("        low_Ea[i]   = low_Ea_DEF[i];")
        self._write("")
        self._write("        rev_A[i]    = rev_A_DEF[i];")
        self._write("        rev_beta[i] = rev_beta_DEF[i];")
        self._write("        rev_Ea[i]   = rev_Ea_DEF[i];")
        self._write("")
        self._write("        troe_a[i]    = troe_a_DEF[i];")
        self._write("        troe_Ts[i]   = troe_Ts_DEF[i];")
        self._write("        troe_Tss[i]  = troe_Tss_DEF[i];")
        self._write("        troe_Tsss[i] = troe_Tsss_DEF[i];")
        self._write("")
        self._write("        sri_a[i] = sri_a_DEF[i];")
        self._write("        sri_b[i] = sri_b_DEF[i];")
        self._write("        sri_c[i] = sri_c_DEF[i];")
        self._write("        sri_d[i] = sri_d_DEF[i];")
        self._write("        sri_e[i] = sri_e_DEF[i];")
        self._write("")
        self._write("        is_PD[i]    = is_PD_DEF[i];")
        self._write("        troe_len[i] = troe_len_DEF[i];")
        self._write("        sri_len[i]  = sri_len_DEF[i];")
        self._write("")
        self._write("        activation_units[i] = activation_units_DEF[i];")
        self._write("        prefactor_units[i]  = prefactor_units_DEF[i];")
        self._write("        phase_units[i]      = phase_units_DEF[i];")
        self._write("")
        self._write("        nTB[i]  = nTB_DEF[i];")
        self._write("        if (nTB[i] != 0) {")
        self._write("           TB[i] = (double *) malloc(sizeof(double) * nTB[i]);")
        self._write("           TBid[i] = (int *) malloc(sizeof(int) * nTB[i]);")
        self._write("           for (int j=0; j<nTB[i]; j++) {")
        self._write("             TB[i][j] = TB_DEF[i][j];")
        self._write("             TBid[i][j] = TBid_DEF[i][j];")
        self._write("           }")
        self._write("        }")
        self._write("    }")
        self._write("}")
        self._write()
        self._write("void SetAllDefaults()")
        self._write("{")
        self._write("    for (int i=0; i<%d; i++) {" % (nReactions))
        self._write("        if (nTB_DEF[i] != 0) {")
        self._write("            nTB_DEF[i] = 0;")
        self._write("            free(TB_DEF[i]);")
        self._write("            free(TBid_DEF[i]);")
        self._write("        }")
        self._write("")
        self._write("        fwd_A_DEF[i]    = fwd_A[i];")
        self._write("        fwd_beta_DEF[i] = fwd_beta[i];")
        self._write("        fwd_Ea_DEF[i]   = fwd_Ea[i];")
        self._write("")
        self._write("        low_A_DEF[i]    = low_A[i];")
        self._write("        low_beta_DEF[i] = low_beta[i];")
        self._write("        low_Ea_DEF[i]   = low_Ea[i];")
        self._write("")
        self._write("        rev_A_DEF[i]    = rev_A[i];")
        self._write("        rev_beta_DEF[i] = rev_beta[i];")
        self._write("        rev_Ea_DEF[i]   = rev_Ea[i];")
        self._write("")
        self._write("        troe_a_DEF[i]    = troe_a[i];")
        self._write("        troe_Ts_DEF[i]   = troe_Ts[i];")
        self._write("        troe_Tss_DEF[i]  = troe_Tss[i];")
        self._write("        troe_Tsss_DEF[i] = troe_Tsss[i];")
        self._write("")
        self._write("        sri_a_DEF[i] = sri_a[i];")
        self._write("        sri_b_DEF[i] = sri_b[i];")
        self._write("        sri_c_DEF[i] = sri_c[i];")
        self._write("        sri_d_DEF[i] = sri_d[i];")
        self._write("        sri_e_DEF[i] = sri_e[i];")
        self._write("")
        self._write("        is_PD_DEF[i]    = is_PD[i];")
        self._write("        troe_len_DEF[i] = troe_len[i];")
        self._write("        sri_len_DEF[i]  = sri_len[i];")
        self._write("")
        self._write("        activation_units_DEF[i] = activation_units[i];")
        self._write("        prefactor_units_DEF[i]  = prefactor_units[i];")
        self._write("        phase_units_DEF[i]      = phase_units[i];")
        self._write("")
        self._write("        nTB_DEF[i]  = nTB[i];")
        self._write("        if (nTB_DEF[i] != 0) {")
        self._write("           TB_DEF[i] = (double *) malloc(sizeof(double) * nTB_DEF[i]);")
        self._write("           TBid_DEF[i] = (int *) malloc(sizeof(int) * nTB_DEF[i]);")
        self._write("           for (int j=0; j<nTB_DEF[i]; j++) {")
        self._write("             TB_DEF[i][j] = TB[i][j];")
        self._write("             TBid_DEF[i][j] = TBid[i][j];")
        self._write("           }")
        self._write("        }")
        self._write("    }")
        self._write("}")
                
        return

    def _renderDocument(self, mechanism, options=None):

        self._setSpecies(mechanism)

        self.reactionIndex = mechanism._sort_reactions()

        self._includes()
        self._declarations(mechanism)
        self._statics(mechanism)
        self._ckinit(mechanism)

        #self._main(mechanism)

        # chemkin wrappers
        self._ckindx(mechanism)
        self._ckxnum(mechanism)
        self._cksnum(mechanism)
        self._cksyme(mechanism)
        self._cksyms(mechanism)
        self._ckrp(mechanism)
        
        self._ckpx(mechanism)
        self._ckpy(mechanism)
        self._vckpy(mechanism)
        self._ckpc(mechanism)
        self._ckrhox(mechanism)
        self._ckrhoy(mechanism)
        self._ckrhoc(mechanism)
        self._ckwt(mechanism)
        self._ckawt(mechanism)
        self._ckmmwy(mechanism)
        self._ckmmwx(mechanism)
        self._ckmmwc(mechanism)
        self._ckytx(mechanism)
        self._vckytx(mechanism)
        self._ckytcp(mechanism)
        self._ckytcr(mechanism)
        self._ckxty(mechanism)
        self._ckxtcp(mechanism)
        self._ckxtcr(mechanism)
        self._ckctx(mechanism)
        self._ckcty(mechanism)
        
        self._ckcpor(mechanism)
        self._ckhort(mechanism)
        self._cksor(mechanism)
        
        self._ckcvml(mechanism)
        self._ckcpml(mechanism)
        self._ckuml(mechanism)
        self._ckhml(mechanism)
        self._ckgml(mechanism)
        self._ckaml(mechanism)
        self._cksml(mechanism)
        
        self._ckcvms(mechanism)
        self._ckcpms(mechanism)
        self._ckums(mechanism)
        self._ckhms(mechanism)
        self._vckhms(mechanism)
        self._ckgms(mechanism)
        self._ckams(mechanism)
        self._cksms(mechanism)

        self._ckcpbl(mechanism)
        self._ckcpbs(mechanism)
        self._ckcvbl(mechanism)
        self._ckcvbs(mechanism)
        
        self._ckhbml(mechanism)
        self._ckhbms(mechanism)
        self._ckubml(mechanism)
        self._ckubms(mechanism)
        self._cksbml(mechanism)
        self._cksbms(mechanism)
        self._ckgbml(mechanism)
        self._ckgbms(mechanism)
        self._ckabml(mechanism)
        self._ckabms(mechanism)

        self._ckwc(mechanism)
        self._ckwyp(mechanism)
        self._ckwxp(mechanism)
        self._ckwyr(mechanism)
        self._vckwyr(mechanism)
        self._ckwxr(mechanism)
        
        self._ckqc(mechanism)
        self._ckkfkr(mechanism)
        self._ckqyp(mechanism)
        self._ckqxp(mechanism)
        self._ckqyr(mechanism)
        self._ckqxr(mechanism)

        self._cknu(mechanism)
        self._ckncf(mechanism)
        
        self._ckabe(mechanism)
        
        self._ckeqc(mechanism)
        self._ckeqyp(mechanism)
        self._ckeqxp(mechanism)
        self._ckeqyr(mechanism)
        self._ckeqxr(mechanism)
        
        # Fuego Functions
        self._productionRate(mechanism)
        self._vproductionRate(mechanism)
        self._DproductionRate(mechanism)
        self._ajac(mechanism)
        self._dthermodT(mechanism)
        self._progressRate(mechanism)
        self._progressRateFR(mechanism)
        self._equilibriumConstants(mechanism)
        self._thermo(mechanism)
        self._molecularWeight(mechanism)
        self._atomicWeight(mechanism)
        self._T_given_ey(mechanism)
        self._T_given_hy(mechanism)
        return


    def _end(self):
        self._timestamp()
        self._rep += self.footer()
        return


    def _includes(self):
        self._rep += [
            '',
            '#include <math.h>',
            '#include <stdio.h>',
            '#include <string.h>',
            '#include <stdlib.h>'
            ]
        return


    def _declarations(self, mechanism):
        self._rep += [
            '',
            '#if defined(BL_FORT_USE_UPPERCASE)',
            '#define CKINDX CKINDX',
            '#define CKINIT CKINIT',
            '#define CKFINALIZE CKFINALIZE',
            '#define CKXNUM CKXNUM',
            '#define CKSYME CKSYME',
            '#define CKSYMS CKSYMS',
            '#define CKRP CKRP',
            '#define CKPX CKPX',
            '#define CKPY CKPY',
            '#define CKPC CKPC',
            '#define CKRHOX CKRHOX',
            '#define CKRHOY CKRHOY',
            '#define CKRHOC CKRHOC',
            '#define CKWT CKWT',
            '#define CKAWT CKAWT',
            '#define CKMMWY CKMMWY',
            '#define CKMMWX CKMMWX',
            '#define CKMMWC CKMMWC',
            '#define CKYTX CKYTX',
            '#define CKYTCP CKYTCP',
            '#define CKYTCR CKYTCR',
            '#define CKXTY CKXTY',
            '#define CKXTCP CKXTCP',
            '#define CKXTCR CKXTCR',
            '#define CKCTX CKCTX',
            '#define CKCTY CKCTY',
            '#define CKCPOR CKCPOR',
            '#define CKHORT CKHORT',
            '#define CKSOR CKSOR',
            '#define CKCVML CKCVML',
            '#define CKCPML CKCPML',
            '#define CKUML CKUML',
            '#define CKHML CKHML',
            '#define CKGML CKGML',
            '#define CKAML CKAML',
            '#define CKSML CKSML',
            '#define CKCVMS CKCVMS',
            '#define CKCPMS CKCPMS',
            '#define CKUMS CKUMS',
            '#define CKHMS CKHMS',
            '#define CKGMS CKGMS',
            '#define CKAMS CKAMS',
            '#define CKSMS CKSMS',
            '#define CKCPBL CKCPBL',
            '#define CKCPBS CKCPBS',
            '#define CKCVBL CKCVBL',
            '#define CKCVBS CKCVBS',
            '#define CKHBML CKHBML',
            '#define CKHBMS CKHBMS',
            '#define CKUBML CKUBML',
            '#define CKUBMS CKUBMS',
            '#define CKSBML CKSBML',
            '#define CKSBMS CKSBMS',
            '#define CKGBML CKGBML',
            '#define CKGBMS CKGBMS',
            '#define CKABML CKABML',
            '#define CKABMS CKABMS',
            '#define CKWC CKWC',
            '#define CKWYP CKWYP',
            '#define CKWXP CKWXP',
            '#define CKWYR CKWYR',
            '#define CKWXR CKWXR',
            '#define CKQC CKQC',
            '#define CKKFKR CKKFKR',
            '#define CKQYP CKQYP',
            '#define CKQXP CKQXP',
            '#define CKQYR CKQYR',
            '#define CKQXR CKQXR',
            '#define CKNU CKNU',
            '#define CKNCF CKNCF',
            '#define CKABE CKABE',
            '#define CKEQC CKEQC',
            '#define CKEQYP CKEQYP',
            '#define CKEQXP CKEQXP',
            '#define CKEQYR CKEQYR',
            '#define CKEQXR CKEQXR',
            '#define DWDOT DWDOT',
            '#define VCKHMS VCKHMS',
            '#define VCKPY VCKPY',
            '#define VCKWYR VCKWYR',
            '#define VCKYTX VCKYTX',
            '#define GET_T_GIVEN_EY GET_T_GIVEN_EY',
            '#define GET_T_GIVEN_HY GET_T_GIVEN_HY',
            '#define GET_REACTION_MAP GET_REACTION_MAP',
            '#elif defined(BL_FORT_USE_LOWERCASE)',
            '#define CKINDX ckindx',
            '#define CKINIT ckinit',
            '#define CKFINALIZE ckfinalize',
            '#define CKXNUM ckxnum',
            '#define CKSYME cksyme',
            '#define CKSYMS cksyms',
            '#define CKRP ckrp',
            '#define CKPX ckpx',
            '#define CKPY ckpy',
            '#define CKPC ckpc',
            '#define CKRHOX ckrhox',
            '#define CKRHOY ckrhoy',
            '#define CKRHOC ckrhoc',
            '#define CKWT ckwt',
            '#define CKAWT ckawt',
            '#define CKMMWY ckmmwy',
            '#define CKMMWX ckmmwx',
            '#define CKMMWC ckmmwc',
            '#define CKYTX ckytx',
            '#define CKYTCP ckytcp',
            '#define CKYTCR ckytcr',
            '#define CKXTY ckxty',
            '#define CKXTCP ckxtcp',
            '#define CKXTCR ckxtcr',
            '#define CKCTX ckctx',
            '#define CKCTY ckcty',
            '#define CKCPOR ckcpor',
            '#define CKHORT ckhort',
            '#define CKSOR cksor',
            '#define CKCVML ckcvml',
            '#define CKCPML ckcpml',
            '#define CKUML ckuml',
            '#define CKHML ckhml',
            '#define CKGML ckgml',
            '#define CKAML ckaml',
            '#define CKSML cksml',
            '#define CKCVMS ckcvms',
            '#define CKCPMS ckcpms',
            '#define CKUMS ckums',
            '#define CKHMS ckhms',
            '#define CKGMS ckgms',
            '#define CKAMS ckams',
            '#define CKSMS cksms',
            '#define CKCPBL ckcpbl',
            '#define CKCPBS ckcpbs',
            '#define CKCVBL ckcvbl',
            '#define CKCVBS ckcvbs',
            '#define CKHBML ckhbml',
            '#define CKHBMS ckhbms',
            '#define CKUBML ckubml',
            '#define CKUBMS ckubms',
            '#define CKSBML cksbml',
            '#define CKSBMS cksbms',
            '#define CKGBML ckgbml',
            '#define CKGBMS ckgbms',
            '#define CKABML ckabml',
            '#define CKABMS ckabms',
            '#define CKWC ckwc',
            '#define CKWYP ckwyp',
            '#define CKWXP ckwxp',
            '#define CKWYR ckwyr',
            '#define CKWXR ckwxr',
            '#define CKQC ckqc',
            '#define CKKFKR ckkfkr',
            '#define CKQYP ckqyp',
            '#define CKQXP ckqxp',
            '#define CKQYR ckqyr',
            '#define CKQXR ckqxr',
            '#define CKNU cknu',
            '#define CKNCF ckncf',
            '#define CKABE ckabe',
            '#define CKEQC ckeqc',
            '#define CKEQYP ckeqyp',
            '#define CKEQXP ckeqxp',
            '#define CKEQYR ckeqyr',
            '#define CKEQXR ckeqxr',
            '#define DWDOT dwdot',
            '#define VCKHMS vckhms',
            '#define VCKPY vckpy',
            '#define VCKWYR vckwyr',
            '#define VCKYTX vckytx',
            '#define GET_T_GIVEN_EY get_t_given_ey',
            '#define GET_T_GIVEN_HY get_t_given_hy',
            '#define GET_REACTION_MAP get_reaction_map',
            '#elif defined(BL_FORT_USE_UNDERSCORE)',
            '#define CKINDX ckindx_',
            '#define CKINIT ckinit_',
            '#define CKFINALIZE ckfinalize_',
            '#define CKXNUM ckxnum_',
            '#define CKSYME cksyme_',
            '#define CKSYMS cksyms_',
            '#define CKRP ckrp_',
            '#define CKPX ckpx_',
            '#define CKPY ckpy_',
            '#define CKPC ckpc_',
            '#define CKRHOX ckrhox_',
            '#define CKRHOY ckrhoy_',
            '#define CKRHOC ckrhoc_',
            '#define CKWT ckwt_',
            '#define CKAWT ckawt_',
            '#define CKMMWY ckmmwy_',
            '#define CKMMWX ckmmwx_',
            '#define CKMMWC ckmmwc_',
            '#define CKYTX ckytx_',
            '#define CKYTCP ckytcp_',
            '#define CKYTCR ckytcr_',
            '#define CKXTY ckxty_',
            '#define CKXTCP ckxtcp_',
            '#define CKXTCR ckxtcr_',
            '#define CKCTX ckctx_',
            '#define CKCTY ckcty_',
            '#define CKCPOR ckcpor_',
            '#define CKHORT ckhort_',
            '#define CKSOR cksor_',
            '#define CKCVML ckcvml_',
            '#define CKCPML ckcpml_',
            '#define CKUML ckuml_',
            '#define CKHML ckhml_',
            '#define CKGML ckgml_',
            '#define CKAML ckaml_',
            '#define CKSML cksml_',
            '#define CKCVMS ckcvms_',
            '#define CKCPMS ckcpms_',
            '#define CKUMS ckums_',
            '#define CKHMS ckhms_',
            '#define CKGMS ckgms_',
            '#define CKAMS ckams_',
            '#define CKSMS cksms_',
            '#define CKCPBL ckcpbl_',
            '#define CKCPBS ckcpbs_',
            '#define CKCVBL ckcvbl_',
            '#define CKCVBS ckcvbs_',
            '#define CKHBML ckhbml_',
            '#define CKHBMS ckhbms_',
            '#define CKUBML ckubml_',
            '#define CKUBMS ckubms_',
            '#define CKSBML cksbml_',
            '#define CKSBMS cksbms_',
            '#define CKGBML ckgbml_',
            '#define CKGBMS ckgbms_',
            '#define CKABML ckabml_',
            '#define CKABMS ckabms_',
            '#define CKWC ckwc_',
            '#define CKWYP ckwyp_',
            '#define CKWXP ckwxp_',
            '#define CKWYR ckwyr_',
            '#define CKWXR ckwxr_',
            '#define CKQC ckqc_',
            '#define CKKFKR ckkfkr_',
            '#define CKQYP ckqyp_',
            '#define CKQXP ckqxp_',
            '#define CKQYR ckqyr_',
            '#define CKQXR ckqxr_',
            '#define CKNU cknu_',
            '#define CKNCF ckncf_',
            '#define CKABE ckabe_',
            '#define CKEQC ckeqc_',
            '#define CKEQYP ckeqyp_',
            '#define CKEQXP ckeqxp_',
            '#define CKEQYR ckeqyr_',
            '#define CKEQXR ckeqxr_',
            '#define DWDOT dwdot_',
            '#define VCKHMS vckhms_',
            '#define VCKPY vckpy_',
            '#define VCKWYR vckwyr_',
            '#define VCKYTX vckytx_',
            '#define GET_T_GIVEN_EY get_t_given_ey_',
            '#define GET_T_GIVEN_HY get_t_given_hy_',
            '#define GET_REACTION_MAP get_reaction_map_',
            '#endif','',
            self.line('function declarations'),
            'void atomicWeight(double * restrict awt);',
            'void molecularWeight(double * restrict wt);',
            'void gibbs(double * restrict species, double * restrict tc);',
            'void helmholtz(double * restrict species, double * restrict tc);',
            'void speciesInternalEnergy(double * restrict species, double * restrict tc);',
            'void speciesEnthalpy(double * restrict species, double * restrict tc);',
            'void speciesEntropy(double * restrict species, double * restrict tc);',
            'void cp_R(double * restrict species, double * restrict tc);',
            'void cv_R(double * restrict species, double * restrict tc);',
            'void equilibriumConstants(double * restrict kc, double * restrict g_RT, double T);',
            'void productionRate(double * restrict wdot, double * restrict sc, double T);',
            'void comp_k_f(double * restrict tc, double invT, double * restrict k_f);',
            'void comp_Kc(double * restrict tc, double invT, double * restrict Kc);',
            'void comp_qfqr(double * restrict q_f, double * restrict q_r, double * restrict sc, double * restrict tc, double invT);',
            'void progressRate(double * restrict qdot, double * restrict speciesConc, double T);',
            'void progressRateFR(double * restrict q_f, double * restrict q_r, double * restrict speciesConc, double T);',
            'void CKINIT'+sym+'();',
            'void CKFINALIZE'+sym+'();',
            'void CKINDX'+sym+'(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit );',
            'void CKXNUM'+sym+'(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline);',
            'void CKSNUM'+sym+'(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict rval, int * kerr, int lenline, int lenkray);',
            'void CKSYME(int * kname, int * lenkname);',
            'void CKSYMS(int * kname, int * lenkname);',
            #'void CKSYMS'+sym+'(char * cckwrk, int * lout, char * kname, int * kerr, int lencck, int lenkname);',
            'void CKRP'+sym+'(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa);',
            'void CKPX'+sym+'(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict P);',
            'void CKPY'+sym+'(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);',
            'void CKPC'+sym+'(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P);',
            'void CKRHOX'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho);',
            'void CKRHOY'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho);',
            'void CKRHOC'+sym+'(double * restrict P, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict rho);',
            'void CKWT'+sym+'(int * iwrk, double * restrict rwrk, double * restrict wt);',
            'void CKAWT'+sym+'(int * iwrk, double * restrict rwrk, double * restrict awt);',
            'void CKMMWY'+sym+'(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm);',
            'void CKMMWX'+sym+'(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm);',
            'void CKMMWC'+sym+'(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm);',
            'void CKYTX'+sym+'(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);',
            'void CKYTCP'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);',
            'void CKYTCR'+sym+'(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);',
            'void CKXTY'+sym+'(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y);',
            'void CKXTCP'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c);',
            'void CKXTCR'+sym+'(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c);',
            'void CKCTX'+sym+'(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x);',
            'void CKCTY'+sym+'(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y);',
            'void CKCPOR'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpor);',
            'void CKHORT'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hort);',
            'void CKSOR'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sor);',
            
            'void CKCVML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml);',
            'void CKCPML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml);',
            'void CKUML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml);',
            'void CKHML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml);',
            'void CKGML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gml);',
            'void CKAML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict aml);',
            'void CKSML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sml);',
            
            'void CKCVMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);',
            'void CKCPMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);',
            'void CKUMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);',
            'void CKHMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);',
            'void CKGMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms);',
            'void CKAMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams);',
            'void CKSMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms);',
            
            'void CKCPBL'+sym+'(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);',
            'void CKCPBS'+sym+'(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);',
            'void CKCVBL'+sym+'(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);',
            'void CKCVBS'+sym+'(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);',
            
            'void CKHBML'+sym+'(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml);',
            'void CKHBMS'+sym+'(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms);',
            'void CKUBML'+sym+'(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml);',
            'void CKUBMS'+sym+'(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms);',
            'void CKSBML'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict sbml);',
            'void CKSBMS'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict sbms);',
            'void CKGBML'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict gbml);',
            'void CKGBMS'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms);',
            'void CKABML'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml);',
            'void CKABMS'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms);',

            
            'void CKWC'+sym+'(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot);',
            'void CKWYP'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot);',
            'void CKWXP'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot);',
            'void CKWYR'+sym+'(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot);',
            'void CKWXR'+sym+'(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot);',

            
            'void CKQC'+sym+'(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot);',
            'void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r);',
            'void CKQYP'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot);',
            'void CKQXP'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot);',
            'void CKQYR'+sym+'(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot);',
            'void CKQXR'+sym+'(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot);',
            
            'void CKNU'+sym+'(int * kdim, int * iwrk, double * restrict rwrk, int * nuki);',
            'void CKNCF'+sym+'(int * mdim, int * iwrk, double * restrict rwrk, int * ncf);',
            
            'void CKABE'+sym+'(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e );',
            'void CKEQC'+sym+'(double * restrict T, double * restrict C , int * iwrk, double * restrict rwrk, double * restrict eqcon );',
            'void CKEQYP'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon);',
            'void CKEQXP'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon);',
            'void CKEQYR'+sym+'(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon);',
            'void CKEQXR'+sym+'(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon);',
            'void DWDOT(double * restrict J, double * restrict sc, double * restrict T, int * consP);',
            'void aJacobian(double * restrict J, double * restrict sc, double T, int consP);',
            'void dcvpRdT(double * restrict species, double * restrict tc);',
            'void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int *ierr);',
            'void GET_T_GIVEN_HY(double * restrict h, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int *ierr);',
            'void GET_REACTION_MAP(int * restrict rmap);',
            self.line('vector version'),
            'void vproductionRate(int npt, double * restrict wdot, double * restrict c, double * restrict T);',
            'void VCKHMS'+sym+'(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);',
            'void VCKPY'+sym+'(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);',
            'void VCKWYR'+sym+'(int * restrict np, double * restrict rho, double * restrict T,',
            '            double * restrict y, int * restrict iwrk, double * restrict rwrk,',
            '            double * restrict wdot);',
            'void VCKYTX'+sym+'(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);',
            'void vcomp_k_f(int npt, double * restrict k_f_s, double * restrict tc, double * restrict invT);',
            'void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc);',
            'void vcomp_Kc(int npt, double * restrict Kc_s, double * restrict g_RT, double * restrict invT);',
            ]
        nReactions = len(mechanism.reaction())
        if nReactions <= 50:
            self._rep += [
                'void vcomp_wdot(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,',
                '                double * restrict k_f_s, double * restrict Kc_s,',
                '                double * restrict tc, double * restrict invT, double * restrict T);',
                ]
        else:
            for i in range(0,nReactions,50):
                self._rep += [
                    'void vcomp_wdot_%d_%d(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,' % (i+1,min(i+50,nReactions)),
                    '                double * restrict k_f_s, double * restrict Kc_s,',
                    '                double * restrict tc, double * restrict invT, double * restrict T);',
                    ]                
        return


    def _main(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('optional test program'))
        self._write('int main()')
        self._write('{')
        self._indent()

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        # declarations
        self._write('int species;')
        self._write('int reaction;')

        self._write('double T;')
        self._write('double q_dot[%d];' % nReactions)
        self._write('double wdot[%d];' % nSpecies)
        self._write('double sc[%d];' % nSpecies)
        self._write('double uml[%d];' % nSpecies)
        self._write('double rckdummy[%d];' % nSpecies)
        self._write('int    ickdummy[%d];' % nSpecies)

        # set the temperature
        self._write()
        self._write('T = 1000.0;')

        # compute ckuml 
        self._write()
        self._write(self.line('compute the internal energy'))
        self._write('CKUML(&T, ickdummy, rckdummy, uml);')
        
        # print
        self._write()
        self._write('for (species = 0; species < %d; ++species) {' % nSpecies)
        self._indent()
        self._write('printf(" e: %5d   %15.7e\\n", species+1, uml[species]);')
        self._outdent()
        self._write('}')


        # compute the gibbs free energy
        # self._write()
        # self._write(self.line('compute the Gibbs free energy'))
        # self._write('gibbs(g_RT, T);')

        # compute the equilibrium constants
        # self._write()
        # self._write(self.line('compute the equilibrium constants'))
        # self._write('equilibriumConstants(kc, g_RT, T);')

        self._write('for (species = 0; species < %d; ++species) {' % nSpecies)
        self._indent()
        self._write('sc[species] = 1.0e6;')
        self._outdent()
        self._write('}')

        # compute the production rates
        self._write()
        self._write(self.line('compute the production rate'))
        self._write('productionRate(wdot, sc, T);')

        # compute the progress rates
        # self._write()
        # self._write(self.line('compute the progress rates'))
        # self._write('progressRate(q_dot, sc, T);')

        # print
        self._write()
        self._write('for (species = 0; species < %d; ++species) {' % nSpecies)
        self._indent()
        self._write('printf("%5d   %15.7e\\n", species+1, wdot[species]);')
        self._outdent()
        self._write('}')

        # print
        # self._write()
        # self._write('for (reaction = 0; reaction < %d; ++reaction) {' % nReactions)
        # self._indent()
        # self._write('printf("%5d | %15.7e\\n", reaction+1, q_dot[reaction]);')
        # self._write('}')
        # self._outdent()

        # done
        self._write()
        self._write('return 0;')

        self._outdent()
        self._write('}')
        return

    def _ckinit(self, mechanism):

        nElement = len(mechanism.element())
        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write(self.line(' Finalizes parameter database'))
        self._write('void CKFINALIZE()')
        self._write('{')
        self._write('  for (int i=0; i<%d; ++i) {' % (nReactions))
        self._write('    free(TB[i]); TB[i] = 0; ')
        self._write('    free(TBid[i]); TBid[i] = 0;')
        self._write('    nTB[i] = 0;')
        self._write()
        self._write('    free(TB_DEF[i]); TB_DEF[i] = 0; ')
        self._write('    free(TBid_DEF[i]); TBid_DEF[i] = 0;')
        self._write('    nTB_DEF[i] = 0;')
        self._write('  }')
        self._write('}')
        self._write()
        self._write(self.line(' Initializes parameter database'))
        self._write('void CKINIT'+sym+'()')
        self._write('{')

        self._indent()

        # build reverse reaction map
        rmap = {}
        for i, reaction in zip(range(nReactions), mechanism.reaction()):
            rmap[reaction.orig_id-1] = i

        for j in range(nReactions):
            reaction = mechanism.reaction()[rmap[j]]
            id = reaction.id - 1

            A, beta, E = reaction.arrhenius
            self._write("// (%d):  %s" % (reaction.orig_id - 1, reaction.equation()))
            self._write("fwd_A[%d]     = %.17g;" % (id,A))
            self._write("fwd_beta[%d]  = %.17g;" % (id,beta))
            self._write("fwd_Ea[%d]    = %.17g;" % (id,E))

            dim = self._phaseSpaceUnits(reaction.reactants)
            thirdBody = reaction.thirdBody
            low = reaction.low
            if not thirdBody:
                uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 3 !PD, !TB
            elif not low:
                uc = self._prefactorUnits(reaction.units["prefactor"], -dim) # Case 2 !PD, TB
            else:
                uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 1 PD, TB
                low_A, low_beta, low_E = low
                self._write("low_A[%d]     = %.17g;" % (id,low_A))
                self._write("low_beta[%d]  = %.17g;" % (id,low_beta))
                self._write("low_Ea[%d]    = %.17g;" % (id,low_E))
                if reaction.troe:
                    troe = reaction.troe
                    ntroe = len(troe)
                    is_troe = True
                    self._write("troe_a[%d]    = %.17g;" % (id,troe[0]))
                    if ntroe>1:
                        self._write("troe_Tsss[%d] = %.17g;" % (id,troe[1]))
                    if ntroe>2:
                        self._write("troe_Ts[%d]   = %.17g;" % (id,troe[2]))
                    if ntroe>3:
                        self._write("troe_Tss[%d]  = %.17g;" % (id,troe[3]))
                    self._write("troe_len[%d]  = %d;" % (id,ntroe))
                if reaction.sri:
                    sri = reaction.sri
                    nsri = len(sri)
                    is_sri = True
                    self._write("sri_a[%d]     = %.17g;" % (id,sri[0]))
                    if nsri>1:
                        self._write("sri_b[%d]     = %.17g;" % (id,sri[1]))
                    if nsri>2:
                        self._write("sri_c[%d]     = %.17g;" % (id,sri[2]))
                    if nsri>3:
                        self._write("sri_d[%d]     = %.17g;" % (id,sri[3]))
                    if nsri>4:
                        self._write("sri_e[%d]     = %.17g;" % (id,sri[4]))
                    self._write("sri_len[%d]   = %d;" % (id,nsri))

            self._write("prefactor_units[%d]  = %.17g;" % (id,uc.value))
            aeuc = self._activationEnergyUnits(reaction.units["activation"])
            self._write("activation_units[%d] = %.17g;" % (id,aeuc / Rc / kelvin))
            self._write("phase_units[%d]      = 1e-%d;" % (id,dim*6))

            if low:
                self._write("is_PD[%d] = 1;" % (id) )
            else:
                self._write("is_PD[%d] = 0;" % (id) )


            if thirdBody:
                efficiencies = reaction.efficiencies
                self._write("nTB[%d] = %d;" % (id, len(efficiencies)))
                self._write("TB[%d] = (double *) malloc(%d * sizeof(double));" % (id, len(efficiencies)))
                self._write("TBid[%d] = (int *) malloc(%d * sizeof(int));" % (id, len(efficiencies)))
                for i, eff in enumerate(efficiencies):
                    symbol, efficiency = eff
                    self._write("TBid[%d][%d] = %.17g; TB[%d][%d] = %.17g; // %s"
                                % (id, i, mechanism.species(symbol).id, id, i, efficiency, symbol ))
            else:
                self._write("nTB[%d] = 0;" % (id))

            self._write()

        self._write("SetAllDefaults();")
        self._outdent()
        self._write("}")
        self._write()
            
        return

    def _thermo(self, mechanism):
        speciesInfo = self._analyzeThermodynamics(mechanism)

        self._gibbs(speciesInfo)
        self._helmholtz(speciesInfo)
        self._cv(speciesInfo)
        self._cp(speciesInfo)
        self._speciesInternalEnergy(speciesInfo)
        self._speciesEnthalpy(speciesInfo)
        self._speciesEntropy(speciesInfo)
        return


    def _dthermodT(self, mechanism):
        speciesInfo = self._analyzeThermodynamics(mechanism)
        self._dcvpdT(speciesInfo)
        return


    def _ckxnum(self, mechanism):
        self._write()
        # self._write()
        # self._write(self.line(' strtok_r: re-entrant (threadsafe) version of strtok, helper function for tokenizing strings '))
        # self._write('char *strtok_r(char *s, const char *delim, char **save_ptr)')
        # self._write('{')
        # self._indent()
        # self._write('char *token;')
        # self._write()
        # self._write('if (s == NULL)')
        # self._indent()
        # self._write('s = *save_ptr;')
        # self._outdent()
        # self._write()
        # self._write('/* Scan leading delimiters.  */')
        # self._write('s += strspn (s, delim);')
        # self._write('if (*s == \'\\0\')')
        # self._write('{')
        # self._indent()
        # self._write('*save_ptr = s;')
        # self._write('return NULL;')
        # self._outdent()
        # self._write('}')
        # self._write()
        # self._write('/* Find the end of the token.  */')
        # self._write('token = s;')
        # self._write('s = strpbrk (token, delim);')
        # self._write('if (s == NULL)')
        # self._indent()
        # self._write('/* This token finishes the string.  */')
        # self._write('*save_ptr = __rawmemchr (token, \'\\0\');')
        # self._outdent()
        # self._write('else')
        # self._write('{')
        # self._indent()
        # self._write('/* Terminate the token and make *SAVE_PTR point past it.  */')
        # self._write('*s = \'\\0\';')
        # self._write('*save_ptr = s + 1;')
        # self._outdent()
        # self._write('}')
        # self._write('return token;')
        # self._outdent()
        # self._write('}')
        # self._write()

        self._write()
        self._write()
        self._write(self.line(' ckxnum... for parsing strings '))
        self._write('void CKXNUM'+sym+'(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline )')
        self._write('{')
        self._indent()
        self._write('int n,i; /*Loop Counters */')
        self._write('char cstr[1000];')
        self._write('char *saveptr;')
        self._write('char *p; /*String Tokens */')
        self._write(self.line(' Strip Comments '))
        self._write('for (i=0; i<lenline; ++i) {')
        self._indent()
        self._write('if (line[i]==\'!\') {')
        self._indent()
        self._write('break;')
        self._outdent()
        self._write('}')
        self._write('cstr[i] = line[i];')
        self._outdent()
        self._write('}')
        self._write('cstr[i] = \'\\0\';')
        self._write()
        self._write('p = strtok_r(cstr," ", &saveptr);')
        self._write('if (!p) {')
        self._indent()
        self._write('*nval = 0;')
        self._write('*kerr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('for (n=0; n<*nexp; ++n) {')
        self._indent()
        self._write('rval[n] = atof(p);')
        self._write('p = strtok_r(NULL, \" \", &saveptr);')
        self._write('if (!p) break;')
        self._outdent()
        self._write('}')
        self._write('*nval = n+1;')
        self._write('if (*nval < *nexp) *kerr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        return

    def _cksnum(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(' cksnum... for parsing strings '))
        self._write('void CKSNUM'+sym+'(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict rval, int * kerr, int lenline, int lenkray)')
        self._write('{')
        self._indent()
        
        self._write(self.line('Not done yet ...'))
        
        # done
        self._outdent()
        self._write('}')
        return

    def _ckrp(self, mechanism):
        self._write()
        self._write()
        self._write(
            self.line(' Returns R, Rc, Patm' ))
        self._write('void CKRP'+sym+'(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa)')
        self._write('{')
        self._indent()
        
        self._write(' *ru  = %g; ' % (R * mole * kelvin / erg))
        self._write(' *ruc = %.20f; ' % (Rc * mole * kelvin / cal))
        self._write(' *pa  = %g; ' % (Patm) )
        
        # done
        self._outdent()
        self._write('}')
        return

    def _cksyme(self, mechanism):

        nElement = len(mechanism.element())
        
        self._write()
        self._write()
        self._write(
            self.line(' Returns the char strings of element names'))
        self._write('void CKSYME'+sym+'(int * kname, int * plenkname )')
        self._write('{')
        self._indent()

        self._write('int i; '+self.line('Loop Counter'))
        self._write('int lenkname = *plenkname;')
        self._write(self.line('clear kname'))
        self._write('for (i=0; i<lenkname*%d; i++) {' % nElement)
        self._indent()
        self._write('kname[i] = \' \';')
        self._outdent()
        self._write('}')
        self._write()
        for element in mechanism.element():
            self._write(self.line(' %s ' % element.symbol))
            ii = 0
            for char in element.symbol:
                self._write('kname[ %d*lenkname + %d ] = \'%s\';' %
                           (element.id, ii, char.capitalize()))
                ii = ii+1
            self._write('kname[ %d*lenkname + %d ] = \' \';' %
                           (element.id, ii))
            self._write()
            
        # done
        self._outdent()
        self._write('}')
        return


    def _cksyms(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        self._write(
            self.line(' Returns the char strings of species names'))
        #self._write('void CKSYMS'+sym+'(char * cckwrk, int * lout, char * kname, int * kerr, int lencck, int lenkname )')
        self._write('void CKSYMS'+sym+'(int * kname, int * plenkname )')
        self._write('{')
        self._indent()
        
        self._write('int i; '+self.line('Loop Counter'))
        self._write('int lenkname = *plenkname;')
        self._write(self.line('clear kname'))
        self._write('for (i=0; i<lenkname*%d; i++) {' % nSpecies)
        self._indent()
        self._write('kname[i] = \' \';')
        self._outdent()
        self._write('}')
        self._write()
        for species in mechanism.species():
            self._write(self.line(' %s ' % species.symbol))
            ii = 0
            for char in species.symbol:
                self._write('kname[ %d*lenkname + %d ] = \'%s\';' %
                           (species.id, ii, char.capitalize()))
                ii = ii+1
            self._write('kname[ %d*lenkname + %d ] = \' \';' %
                           (species.id, ii))
            self._write()

        # done
        self._outdent()
        self._write('}')
        return


    def _ckindx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('A few mechanism parameters'))
        self._write('void CKINDX'+sym+'(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)')
        self._write('{')
        self._indent()
        self._write('*mm = %d;' % len(mechanism.element()))
        self._write('*kk = %d;' % len(mechanism.species()))
        self._write('*ii = %d;' % len(mechanism.reaction()))
        self._write('*nfit = -1; ' + self.line(
            'Why do you need this anyway ? '))
        
        # done
        self._outdent()
        self._write('}')
        return
        
        
    def _ckpx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute P = rhoRT/W(x)'))
        self._write('void CKPX'+sym+'(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict P)')
        self._write('{')
        self._indent()

        self._write('double XW = 0;'+
                    self.line(' To hold mean molecular wt'))
        
        # molecular weights of all species
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write(
            '*P = *rho * %g * (*T) / XW; ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        return

    def _ckpy(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute P = rhoRT/W(y)'))
        self._write('void CKPY'+sym+'(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)')
        self._write('{')
        self._indent()

        self._write('double YOW = 0;'+self.line(' for computing mean MW'))
        
        # molecular weights of all species
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))

        self.line('YOW holds the reciprocal of the mean molecular wt')
        self._write(
            '*P = *rho * %g * (*T) * YOW; ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 

    def _vckpy(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute P = rhoRT/W(y)'))
        self._write('void VCKPY'+sym+'(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)')
        self._write('{')
        self._indent()

        species = self.species
        nSpec = len(species)
        self._write('double YOW[*np];')
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('YOW[i] = 0.0;')
        self._outdent()
        self._write('}')        
        self._write('')
        self._write('for (int n=0; n<%d; n++) {' % (nSpec))
        self._indent()
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('YOW[i] += y[n*(*np)+i] * imw[n];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write('')

        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write(
            'P[i] = rho[i] * %g * T[i] * YOW[i]; ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        self._outdent()
        self._write('}')
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckpc(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute P = rhoRT/W(c)'))
        self._write('void CKPC'+sym+'(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P)')
        
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write(self.line('See Eq 5 in CK Manual'))
        self._write('double W = 0;')
        self._write('double sumC = 0;')
        
        # molecular weights of all species
        for species in self.species:
            self._write('W += c[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        nSpecies = len(mechanism.species())
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('sumC += c[id];')
        self._outdent()
        self._write('}')

        self.line('W/sumC holds the mean molecular wt')
        self._write(
            '*P = *rho * %g * (*T) * sumC / W; ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 

    def _ckrhox(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute rho = PW(x)/RT'))
        self._write('void CKRHOX'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho)')
        self._write('{')
        self._indent()

        self._write('double XW = 0;'+
                    self.line(' To hold mean molecular wt'))
        
        # molecular weights of all species
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write(
            '*rho = *P * XW / (%g * (*T)); ' % (R*kelvin*mole/erg)
            + self.line('rho = P*W/(R*T)'))
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        return

    def _ckrhoy(self, mechanism):
        species = self.species
        nSpec = len(species)
        self._write()
        self._write()
        self._write(self.line('Compute rho = P*W(y)/RT'))
        self._write('void CKRHOY'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)')
        self._write('{')
        self._indent()
        self._write('double YOW = 0;')
        self._write('double tmp[%d];' % (nSpec))
        self._write('')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('tmp[i] = y[i]*imw[i];')
        self._outdent()
        self._write('}')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('YOW += tmp[i];')
        self._outdent()
        self._write('}')
        self._write('')
        self._write('*rho = *P / (%g * (*T) * YOW);' % (R * mole * kelvin / erg) + self.line('rho = P*W/(R*T)'))
        self._write('return;')
        self._outdent()
        self._write('}')
        return 
 
    def _ckrhoc(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute rho = P*W(c)/(R*T)'))
        self._write('void CKRHOC'+sym+'(double * restrict P, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict rho)')
        
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write(self.line('See Eq 5 in CK Manual'))
        self._write('double W = 0;')
        self._write('double sumC = 0;')
        
        # molecular weights of all species
        for species in self.species:
            self._write('W += c[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('sumC += c[id];')
        self._outdent()
        self._write('}')

        self.line('W/sumC holds the mean molecular wt')
        self._write(
            '*rho = *P * W / (sumC * (*T) * %g); ' % (R*kelvin*mole/erg)
            + self.line('rho = PW/(R*T)'))
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 

    def _ckwt(self, mechanism):

        self._write()
        self._write()
        self._write(self.line('get molecular weight for all species'))
        self._write('void CKWT'+sym+'(int * iwrk, double * restrict rwrk, double * restrict wt)')
        self._write('{')
        self._indent()

        # call molecularWeight
        self._write('molecularWeight(wt);')
        
        self._outdent()

        self._write('}')

        return
      
    def _ckawt(self, mechanism):

        self._write()
        self._write()
        self._write(self.line('get atomic weight for all elements'))
        self._write('void CKAWT'+sym+'(int * iwrk, double * restrict rwrk, double * restrict awt)')
        self._write('{')
        self._indent()

        # call atomicWeight
        self._write('atomicWeight(awt);')
        
        self._outdent()

        self._write('}')

        return
      
    def _ckcvml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get specific heat at constant volume as a function '))
        self._write(self.line('of T for all species (molar units)'))
        self._write('void CKCVML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('cv_R(cvml, tc);')
        
        # convert cv/R to cv
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('cvml[id] *= %g;' % (R*kelvin*mole/erg) )
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return
       
    def _ckcpml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get specific heat at constant pressure as a '))
        self._write(self.line('function of T for all species (molar units)'))
        self._write('void CKCPML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('cp_R(cpml, tc);')
        
        # convert cp/R to cp
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('cpml[id] *= %g;' % (R*kelvin*mole/erg) )
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return
     
    def _ckuml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get internal energy as a function '))
        self._write(self.line('of T for all species (molar units)'))
        self._write('void CKUML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesInternalEnergy(uml, tc);')
        
        # convert e/RT to e with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('uml[id] *= RT;')
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return
      
    def _ckhml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get enthalpy as a function '))
        self._write(self.line('of T for all species (molar units)'))
        self._write('void CKHML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesEnthalpy(hml, tc);')
        
        # convert h/RT to h with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('hml[id] *= RT;')
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return
    
    def _ckgml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get standard-state Gibbs energy as a function '))
        self._write(self.line('of T for all species (molar units)'))
        self._write('void CKGML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('gibbs(gml, tc);')
        
        # convert g/RT to g with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('gml[id] *= RT;')
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return
    
    def _ckaml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get standard-state Helmholtz free energy as a '))
        self._write(self.line('function of T for all species (molar units)'))
        self._write('void CKAML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict aml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('helmholtz(aml, tc);')
        
        # convert A/RT to A with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('aml[id] *= RT;')
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return
   
    def _cksml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the standard-state entropies in molar units'))
        self._write('void CKSML'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('speciesEntropy(sml, tc);')
        
        # convert s/R to s
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('sml[id] *= %g;' % (R*kelvin*mole/erg) )
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return
 
    def _ckums(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns internal energy in mass units (Eq 30.)'))
        self._write('void CKUMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesInternalEnergy(ums, tc);')
        
        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('ums[i] *= RT*imw[i];')
        self._outdent()
        self._write('}')
        self._outdent()

        self._write('}')

        return
 
    def _ckhms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns enthalpy in mass units (Eq 27.)'))
        self._write('void CKHMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesEnthalpy(hms, tc);')
        
        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('hms[i] *= RT*imw[i];')
        self._outdent()
        self._write('}')
        self._outdent()

        self._write('}')

        return

    def _vckhms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns enthalpy in mass units (Eq 27.)'))
        self._write('void VCKHMS'+sym+'(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)')
        self._write('{')
        self._indent()

        species = self.species
        nSpec = len(species)

        self._write('double tc[5], h[%d];' % nSpec)

        self._write()

        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('tc[0] = 0.0;')
        self._write('tc[1] = T[i];')
        self._write('tc[2] = T[i]*T[i];')
        self._write('tc[3] = T[i]*T[i]*T[i];')
        self._write('tc[4] = T[i]*T[i]*T[i]*T[i];')

        self._write()

        self._write('speciesEnthalpy(h, tc);')

        self._write()

        for ispec in range(nSpec):
            self._write('hms[%d*(*np)+i] = h[%d];' % (ispec, ispec))
        self._outdent()
        self._write('}')

        self._write()
        
        self._write('for (int n=0; n<%d; n++) {' % (nSpec))
        self._indent()
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('hms[n*(*np)+i] *= %g * T[i] * imw[n];' % (R*kelvin*mole/erg))
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._outdent()
        self._write('}')

        return

    def _ckams(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns helmholtz in mass units (Eq 32.)'))
        self._write('void CKAMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('helmholtz(ams, tc);')
        
        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('ams[i] *= RT*imw[i];')
        self._outdent()
        self._write('}')
        self._outdent()

        self._write('}')

        return

    def _ckgms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns gibbs in mass units (Eq 31.)'))
        self._write('void CKGMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('gibbs(gms, tc);')
        
        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('gms[i] *= RT*imw[i];')
        self._outdent()
        self._write('}')
        self._outdent()

        self._write('}')

        return


    def _ckcvms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the specific heats at constant volume'))
        self._write(self.line('in mass units (Eq. 29)'))
        self._write('void CKCVMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('cv_R(cvms, tc);')

        # convert cv/R to cv with mass units
        self._write(self.line('multiply by R/molecularweight'))
        for species in self.species:
            ROW = (R*kelvin*mole/erg) / species.weight
            self._write('cvms[%d] *= %20.15e; ' % (
                species.id, ROW) + self.line('%s' % species.symbol))

       
        self._outdent()

        self._write('}')

        return

    def _ckcpms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the specific heats at constant pressure'))
        self._write(self.line('in mass units (Eq. 26)'))
        self._write('void CKCPMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpms)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('cp_R(cpms, tc);')
        

        # convert cp/R to cp with mass units
        self._write(self.line('multiply by R/molecularweight'))
        for species in self.species:
            ROW = (R*kelvin*mole/erg) / species.weight
            self._write('cpms[%d] *= %20.15e; ' % (
                species.id, ROW) + self.line('%s' % species.symbol))

       
        self._outdent()

        self._write('}')

        return

    def _cksms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the entropies in mass units (Eq 28.)'))
        self._write('void CKSMS'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('speciesEntropy(sms, tc);')
        

        # convert s/R to s with mass units
        self._write(self.line('multiply by R/molecularweight'))
        for species in self.species:
            ROW = (R*kelvin*mole/erg) / species.weight
            self._write('sms[%d] *= %20.15e; ' % (
                species.id, ROW) + self.line('%s' % species.symbol))

       
        self._outdent()

        self._write('}')

        return
    
    def _ckcpbl(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the mean specific heat at CP (Eq. 33)'))
        self._write('void CKCPBL'+sym+'(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double cpor[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        # call routine
        self._write('cp_R(cpor, tc);')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*cpor[id];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*cpbl = result * %g;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return
 
    def _ckcpbs(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the mean specific heat at CP (Eq. 34)'))
        self._write('void CKCPBS'+sym+'(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs)')
        self._write('{')
        self._indent()

        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double cpor[%d], tresult[%d]; ' % (self.nSpecies,self.nSpecies) + self.line(' temporary storage'))
        
        # call routine
        self._write('cp_R(cpor, tc);')
        

        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('tresult[i] = cpor[i]*y[i]*imw[i];')
        self._outdent()
        self._write('')
        self._write('}')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('result += tresult[i];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*cpbs = result * %g;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return
   
    def _ckcvbl(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the mean specific heat at CV (Eq. 35)'))
        self._write('void CKCVBL'+sym+'(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cvbl)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double cvor[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        # call routine
        self._write('cv_R(cvor, tc);')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*cvor[id];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*cvbl = result * %g;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return

    def _ckcvbs(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the mean specific heat at CV (Eq. 36)'))
        self._write('void CKCVBS'+sym+'(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cvbs)')
        self._write('{')
        self._indent()

        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double cvor[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        # call routine
        self._write('cv_R(cvor, tc);')
        
        # do dot product
        self._write(self.line('multiply by y/molecularweight'))
        for species in self.species:
            self._write('result += cvor[%d]*y[%d]*imw[%d]; ' % (
                species.id, species.id, species.id) + self.line('%s' % species.symbol))

        self._write()
        self._write('*cvbs = result * %g;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return
    
    def _ckhbml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the mean enthalpy of the mixture in molar units'))
        self._write('void CKHBML'+sym+'(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double hml[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesEnthalpy(hml, tc);')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*hml[id];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*hbml = result * RT;')
        
        self._outdent()

        self._write('}')

        return
 
 
    def _ckhbms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns mean enthalpy of mixture in mass units'))
        self._write('void CKHBMS'+sym+'(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms)')
        self._write('{')
        self._indent()

        self._write('double result = 0;')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double hml[%d], tmp[%d]; ' % (self.nSpecies,self.nSpecies) + self.line(' temporary storage'))
        
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesEnthalpy(hml, tc);')

        self._write('int id;')
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('tmp[id] = y[id]*hml[id]*imw[id];')
        self._outdent()
        self._write('}')
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += tmp[id];')
        self._outdent()
        self._write('}')

        self._write()
        # finally, multiply by RT
        self._write('*hbms = result * RT;')
        
        self._outdent()

        self._write('}')
        
        return
    
    def _ckubml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get mean internal energy in molar units'))
        self._write('void CKUBML'+sym+'(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double uml[%d]; ' % self.nSpecies + self.line(' temporary energy array'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesInternalEnergy(uml, tc);')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*uml[id];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*ubml = result * RT;')
        
        self._outdent()

        self._write('}')

        return
 
    def _ckubms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get mean internal energy in mass units'))
        self._write('void CKUBMS'+sym+'(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms)')
        self._write('{')
        self._indent()

        self._write('double result = 0;')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double ums[%d]; ' % self.nSpecies + self.line(' temporary energy array'))
        
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesInternalEnergy(ums, tc);')

        # convert e/RT to e with mass units
        self._write(self.line('perform dot product + scaling by wt'))
        for species in self.species:
            self._write('result += y[%d]*ums[%d]*imw[%d]; ' % (
                species.id, species.id, species.id)
                        + self.line('%s' % species.symbol))

        
        self._write()
        # finally, multiply by RT
        self._write('*ubms = result * RT;')
        
        self._outdent()

        self._write('}')
        
        return
 
    def _cksbml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get mixture entropy in molar units'))
        self._write('void CKSBML'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict sbml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double sor[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        
        # call routine
        self._write('speciesEntropy(sor, tc);')
        
        # Equation 42
        self._write()
        self._write(self.line('Compute Eq 42'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*(sor[id]-log((x[id]+%g))-logPratio);' %
                    smallnum )
        self._outdent()
        self._write('}')

        self._write()
        
        self._write('*sbml = result * %g;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return

    def _cksbms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get mixture entropy in mass units'))
        self._write('void CKSBMS'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict sbms)')
        self._write('{')
        self._indent()

        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double sor[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        self._write(
            'double x[%d]; ' % self.nSpecies + self.line(' need a ytx conversion'))

        self._write('double YOW = 0; '+self.line('See Eq 4, 6 in CK Manual'))
        
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))
 
        # now to ytx
        self._write(self.line('Now compute y to x conversion'))
        for species in self.species:
            self._write('x[%d] = y[%d]/(%f*YOW); ' % (
                species.id, species.id, species.weight) )
            
        # call routine
        self._write('speciesEntropy(sor, tc);')
        
        # Equation 42 and 43
        self._write(self.line('Perform computation in Eq 42 and 43'))
        for species in self.species:
            self._write('result += x[%d]*(sor[%d]-log((x[%d]+%g))-logPratio);' %
                        (species.id, species.id, species.id, smallnum) )

        self._write(self.line('Scale by R/W'))
        self._write('*sbms = result * %g * YOW;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return

    def _ckgbml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns mean gibbs free energy in molar units'))
        self._write('void CKGBML'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict gbml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        self._write(
            'double gort[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        # call routine
        self._write(self.line('Compute g/RT'))
        self._write('gibbs(gort, tc);')
        
        # Equation 44
        self._write()
        self._write(self.line('Compute Eq 44'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*(gort[id]+log((x[id]+%g))+logPratio);' %
                    smallnum )
        self._outdent()
        self._write('}')

        self._write()
        
        self._write('*gbml = result * RT;')
        
        self._outdent()

        self._write('}')

        return


    def _ckgbms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns mixture gibbs free energy in mass units'))
        self._write('void CKGBMS'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms)')
        self._write('{')
        self._indent()

        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        self._write(
            'double gort[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        self._write(
            'double x[%d]; ' % self.nSpecies + self.line(' need a ytx conversion'))

        self._write(
            'double YOW = 0; '
            + self.line('To hold 1/molecularweight'))
        
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))
 
        # now to ytx
        self._write(self.line('Now compute y to x conversion'))
        for species in self.species:
            self._write('x[%d] = y[%d]/(%f*YOW); ' % (
                species.id, species.id, species.weight) )
            
        # call routine
        self._write('gibbs(gort, tc);')
        
        # Equation 42 and 43
        self._write(self.line('Perform computation in Eq 44'))
        for species in self.species:
            self._write('result += x[%d]*(gort[%d]+log((x[%d]+%g))+logPratio);' %
                        (species.id, species.id, species.id, smallnum) )

        self._write(self.line('Scale by RT/W'))
        self._write('*gbms = result * RT * YOW;')
        
        self._outdent()

        self._write('}')

        return
    

    def _ckabml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns mean helmholtz free energy in molar units'))
        self._write('void CKABML'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        self._write(
            'double aort[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        # call routine
        self._write(self.line('Compute g/RT'))
        self._write('helmholtz(aort, tc);')
        
        # Equation 44
        self._write()
        self._write(self.line('Compute Eq 44'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*(aort[id]+log((x[id]+%g))+logPratio);' %
                    smallnum )
        self._outdent()
        self._write('}')

        self._write()
        
        self._write('*abml = result * RT;')
        
        self._outdent()

        self._write('}')

        return
    

    def _ckabms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns mixture helmholtz free energy in mass units'))
        self._write('void CKABMS'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms)')
        self._write('{')
        self._indent()

        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        self._write(
            'double aort[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        self._write(
            'double x[%d]; ' % self.nSpecies + self.line(' need a ytx conversion'))

        self._write(
            'double YOW = 0; '
            + self.line('To hold 1/molecularweight'))
        
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))
 
        # now to ytx
        self._write(self.line('Now compute y to x conversion'))
        for species in self.species:
            self._write('x[%d] = y[%d]/(%f*YOW); ' % (
                species.id, species.id, species.weight) )
            
        # call routine
        self._write('helmholtz(aort, tc);')
        
        # Equation 42 and 43
        self._write(self.line('Perform computation in Eq 44'))
        for species in self.species:
            self._write('result += x[%d]*(aort[%d]+log((x[%d]+%g))+logPratio);' %
                        (species.id, species.id, species.id, smallnum) )

        self._write(self.line('Scale by RT/W'))
        self._write('*abms = result * RT * YOW;')
        
        self._outdent()

        self._write('}')

        return
    

    def _ckwc(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('compute the production rate for each species'))
        self._write('void CKWC'+sym+'(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        # convert C to SI units
        self._write()
        self._write(self.line('convert to SI'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e6;')
        self._outdent()
        self._write('}')
        
        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('productionRate(wdot, C, *T);')

        # convert C and wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e-6;')
        self._write('wdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return

    def _ckwyp(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the molar production rate of species'))
        self._write(self.line('Given P, T, and mass fractions'))
        self._write('void CKWYP'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % self.nSpecies + self.line('temporary storage'))
        self._write('double YOW = 0; ')
        self._write('double PWORT; ')
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))
 
        self._write(self.line('PW/RT (see Eq. 7)'))
        self._write('PWORT = (*P)/(YOW * %g * (*T)); ' % (R*kelvin*mole/erg) )
        
        self._write(self.line('multiply by 1e6 so c goes to SI'))
        self._write('PWORT *= 1e6; ')

        # now compute conversion
        self._write(self.line('Now compute conversion (and go to SI)'))
        for species in self.species:
            self._write('c[%d] = PWORT * y[%d]*imw[%d]; ' % (
                species.id, species.id, species.id) )

        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('productionRate(wdot, c, *T);')

        # convert wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('wdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckwxp(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the molar production rate of species'))
        self._write(self.line('Given P, T, and mole fractions'))
        self._write('void CKWXP'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % self.nSpecies + self.line('temporary storage'))
        
        self._write('double PORT = 1e6 * (*P)/(%g * (*T)); ' % (R*kelvin*mole/erg) +
                    self.line('1e6 * P/RT so c goes to SI units'))
        
        # now compute conversion
        self._write()
        self._write(self.line('Compute conversion, see Eq 10'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('c[id] = x[id]*PORT;')
        self._outdent()
        self._write('}')
        
        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('productionRate(wdot, c, *T);')

        # convert wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('wdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckwyr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the molar production rate of species'))
        self._write(self.line('Given rho, T, and mass fractions'))
        self._write('void CKWYR'+sym+'(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % self.nSpecies + self.line('temporary storage'))

        # now compute conversion
        self._write(self.line('See Eq 8 with an extra 1e6 so c goes to SI'))
        for species in self.species:
            self._write('c[%d] = 1e6 * (*rho) * y[%d]*imw[%d]; ' % (
                species.id, species.id, species.id) )
            
        # call productionRate
        self._write()
        self._write(self.line('call productionRate'))
        self._write('productionRate(wdot, c, *T);')

        # convert wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('wdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _vckwyr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the molar production rate of species'))
        self._write(self.line('Given rho, T, and mass fractions'))
        self._write('void VCKWYR'+sym+'(int * restrict np, double * restrict rho, double * restrict T,')
        self._write('	    double * restrict y, int * restrict iwrk, double * restrict rwrk,')
        self._write('	    double * restrict wdot)')
        self._write('{')
        self._indent()

        self._write('double c[%d*(*np)]; ' % self.nSpecies + self.line('temporary storage'))

        # now compute conversion
        self._write(self.line('See Eq 8 with an extra 1e6 so c goes to SI'))
        self._write('for (int n=0; n<%d; n++) {' % self.nSpecies)
        self._indent()
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        # call productionRate
        self._write()
        self._write(self.line('call productionRate'))
        self._write('vproductionRate(*np, wdot, c, T);')

        # convert wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (int i=0; i<%d*(*np); i++) {' % self.nSpecies)
        self._indent()
        self._write('wdot[i] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckwxr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the molar production rate of species'))
        self._write(self.line('Given rho, T, and mole fractions'))
        self._write('void CKWXR'+sym+'(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % self.nSpecies + self.line('temporary storage'))
        
        self._write('double XW = 0; '+self.line('See Eq 4, 11 in CK Manual'))
        self._write('double ROW; ')
        
        # compute mean molecular weight first (eq 3)
        self._write(self.line('Compute mean molecular wt first'))
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        # now compute conversion
        self._write(self.line('Extra 1e6 factor to take c to SI'))
        self._write('ROW = 1e6*(*rho) / XW;')
        self._write()
        self._write(self.line('Compute conversion, see Eq 11'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('c[id] = x[id]*ROW;')
        self._outdent()
        self._write('}')
        
        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('productionRate(wdot, c, *T);')

        # convert wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('wdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _cknu(self, mechanism):

        nSpecies  = len(mechanism.species())
        nReaction = len(mechanism.reaction())

        self._write()
        self._write()
        self._write(self.line('Returns the stoichiometric coefficients'))
        self._write(self.line('of the reaction mechanism. (Eq 50)'))
        self._write('void CKNU'+sym+'(int * kdim, int * iwrk, double * restrict rwrk, int * nuki)')
        self._write('{')
        self._indent()

 
        self._write('int id; ' + self.line('loop counter'))
        self._write('int kd = (*kdim); ')
        self._write(self.line('Zero nuki'))
        self._write('for (id = 0; id < %d * kd; ++ id) {' % (nSpecies) )
        self._indent()
        self._write(' nuki[id] = 0; ')
        self._outdent()
        self._write('}')
        
        for reaction in mechanism.reaction():

            self._write()
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

            for symbol, coefficient in reaction.reactants:
                self._write(
                    "nuki[ %d * kd + %d ] += -%d ;"
                    % (mechanism.species(symbol).id, reaction.id-1, coefficient))

            for symbol, coefficient in reaction.products:
                self._write(
                    "nuki[ %d * kd + %d ] += +%d ;"
                    % (mechanism.species(symbol).id, reaction.id-1, coefficient))
       
        # done
        self._outdent()
        self._write('}')

        return


    def _ckncf(self, mechanism):

        nSpecies  = len(mechanism.species())
        nElement  = len(mechanism.element())

        self._write()
        self._write()
        self._write(self.line('Returns the elemental composition '))
        self._write(self.line('of the speciesi (mdim is num of elements)'))
        self._write('void CKNCF'+sym+'(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)')
        self._write('{')
        self._indent()

 
        self._write('int id; ' + self.line('loop counter'))
        self._write('int kd = (*mdim); ')
        self._write(self.line('Zero ncf'))
        self._write('for (id = 0; id < kd * %d; ++ id) {' % (self.nSpecies) )
        self._indent()
        self._write(' ncf[id] = 0; ')
        self._outdent()
        self._write('}')
        
        self._write()
        for species in mechanism.species():
           self._write(self.line('%s' % species.symbol))
           for elem, coef in species.composition:
               self._write('ncf[ %d * kd + %d ] = %d; ' % (
                   species.id, mechanism.element(elem).id, coef) +
                       self.line('%s' % elem) )
                           
           self._write()
                            
        # done
        self._outdent()

        self._write('}')

        return


    def _ckabe(self, mechanism):

        nElement  = len(mechanism.element())

        self._write()
        self._write()
        self._write(self.line('Returns the arrehenius coefficients '))
        self._write(self.line('for all reactions'))
        self._write('void CKABE'+sym+'(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)')
        self._write('{')
        self._indent()

        self._write('for (int i=0; i<%d; ++i) {' % len(mechanism.reaction()) )
        self._indent()
        self._write("a[i] = fwd_A[i];")
        self._write("b[i] = fwd_beta[i];")
        self._write("e[i] = fwd_Ea[i];")
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return
                            
    
    def _ckmmwy(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('given y[species]: mass fractions'))
        self._write(self.line('returns mean molecular weight (gm/mole)'))
        self._write('void CKMMWY'+sym+'(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm)')
        self._write('{')
        self._indent()
        species = self.species
        nSpec = len(species)
        self._write('double YOW = 0;')
        self._write('double tmp[%d];' % (nSpec))
        self._write('')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('tmp[i] = y[i]*imw[i];')
        self._outdent()
        self._write('}')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('YOW += tmp[i];')
        self._outdent()
        self._write('}')
        self._write('')
        self._write('*wtm = 1.0 / YOW;')
        self._write('return;')
        self._outdent()
        self._write('}')
        return 
 
    def _ckmmwx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('given x[species]: mole fractions'))
        self._write(self.line('returns mean molecular weight (gm/mole)'))
        self._write('void CKMMWX'+sym+'(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm)')
        self._write('{')
        self._indent()

        self._write('double XW = 0;'+self.line(' see Eq 4 in CK Manual'))
        
        # molecular weights of all species
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write('*wtm = XW;')
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckmmwc(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('given c[species]: molar concentration'))
        self._write(self.line('returns mean molecular weight (gm/mole)'))
        self._write('void CKMMWC'+sym+'(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write(self.line('See Eq 5 in CK Manual'))
        self._write('double W = 0;')
        self._write('double sumC = 0;')
        
        # molecular weights of all species
        for species in self.species:
            self._write('W += c[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('sumC += c[id];')
        self._outdent()
        self._write('}')

        self._write(self.line(' CK provides no guard against divison by zero'))
        self._write('*wtm = W/sumC;')
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckytx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert y[species] (mass fracs) to x[species] (mole fracs)'))
        self._write('void CKYTX'+sym+'(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)')
        self._write('{')
        self._indent()

        species = self.species
        nSpec = len(species)
        self._write('double YOW = 0;')
        self._write('double tmp[%d];' % (nSpec))
        self._write('')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('tmp[i] = y[i]*imw[i];')
        self._outdent()
        self._write('}')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('YOW += tmp[i];')
        self._outdent()
        self._write('}')
        self._write('')
        self._write('double YOWINV = 1.0/YOW;')
        self._write('')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('x[i] = y[i]*imw[i]*YOWINV;')
        self._outdent()
        self._write('}')
        self._write('return;')
        self._outdent()
        self._write('}')
        return 
 
    def _vckytx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs)'))
        self._write('void VCKYTX'+sym+'(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)')
        self._write('{')
        self._indent()

        species = self.species
        nSpec = len(species)
        self._write('double YOW[*np];')
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('YOW[i] = 0.0;')
        self._outdent()
        self._write('}')        
        self._write('')
        self._write('for (int n=0; n<%d; n++) {' % (nSpec))
        self._indent()
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('x[n*(*np)+i] = y[n*(*np)+i] * imw[n];')
        self._write('YOW[i] += x[n*(*np)+i];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write('')

        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('YOW[i] = 1.0/YOW[i];')
        self._outdent()
        self._write('}')

        self._write('')
        
        self._write('for (int n=0; n<%d; n++) {' % (nSpec))
        self._indent()
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('x[n*(*np)+i] *=  YOW[i];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        return 
 
    def _ckytcp(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert y[species] (mass fracs) to c[species] (molar conc)'))
        self._write('void CKYTCP'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)')
        self._write('{')
        self._indent()


        species = self.species
        nSpec = len(species)
        self._write('double YOW = 0;')
        self._write('double PWORT;')
        self._write('')
        self._write(self.line('Compute inverse of mean molecular wt first'))
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('c[i] = y[i]*imw[i];')
        self._outdent()
        self._write('}')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('YOW += c[i];')
        self._outdent()
        self._write('}')
        self._write('')
        self._write(self.line('PW/RT (see Eq. 7)'))
        self._write('PWORT = (*P)/(YOW * %g * (*T)); ' % (R*kelvin*mole/erg) )

        # now compute conversion
        self._write(self.line('Now compute conversion'))
        self._write('')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('c[i] = PWORT * y[i] * imw[i];')
        self._outdent()
        self._write('}')
        self._write('return;')
        self._outdent()
        self._write('}')
        return 
 
    def _ckytcr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert y[species] (mass fracs) to c[species] (molar conc)'))
        self._write('void CKYTCR'+sym+'(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)')
        self._write('{')
        self._indent()
        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('c[i] = (*rho)  * y[i] * imw[i];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        return 
 
    def _ckxty(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert x[species] (mole fracs) to y[species] (mass fracs)'))
        self._write('void CKXTY'+sym+'(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y)')
        self._write('{')
        self._indent()

        self._write('double XW = 0; '+self.line('See Eq 4, 9 in CK Manual'))
        
        # compute mean molecular weight first (eq 3)
        self._write(self.line('Compute mean molecular wt first'))
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))
 
        # now compute conversion
        self._write(self.line('Now compute conversion'))
        self._write('double XWinv = 1.0/XW;')
        for species in self.species:
            self._write('y[%d] = x[%d]*%f*XWinv; ' % (
                species.id, species.id, species.weight) )

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckxtcp(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert x[species] (mole fracs) to c[species] (molar conc)'))
        self._write('void CKXTCP'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double PORT = (*P)/(%g * (*T)); ' % (R*kelvin*mole/erg) +
                    self.line('P/RT'))
        # now compute conversion
        self._write()
        self._write(self.line('Compute conversion, see Eq 10'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('c[id] = x[id]*PORT;')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckxtcr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert x[species] (mole fracs) to c[species] (molar conc)'))
        self._write('void CKXTCR'+sym+'(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double XW = 0; '+self.line('See Eq 4, 11 in CK Manual'))
        self._write('double ROW; ')
        
        # compute mean molecular weight first (eq 3)
        self._write(self.line('Compute mean molecular wt first'))
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        # now compute conversion
        self._write('ROW = (*rho) / XW;')
        self._write()
        self._write(self.line('Compute conversion, see Eq 11'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('c[id] = x[id]*ROW;')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckctx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert c[species] (molar conc) to x[species] (mole fracs)'))
        self._write('void CKCTX'+sym+'(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double sumC = 0; ')

        self._write()
        self._write(self.line('compute sum of c '))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('sumC += c[id];')
        self._outdent()
        self._write('}')

        # now compute conversion
        self._write()
        self._write(self.line(' See Eq 13 '))
        self._write('double sumCinv = 1.0/sumC;')
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('x[id] = c[id]*sumCinv;')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckcty(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert c[species] (molar conc) to y[species] (mass fracs)'))
        self._write('void CKCTY'+sym+'(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y)')
        self._write('{')
        self._indent()

        self._write('double CW = 0; '+self.line('See Eq 12 in CK Manual'))
        
        # compute denominator in eq 12
        self._write(self.line('compute denominator in eq 12 first'))
        for species in self.species:
            self._write('CW += c[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        # now compute conversion
        self._write(self.line('Now compute conversion'))
        self._write('double CWinv = 1.0/CW;')
        for species in self.species:
            self._write('y[%d] = c[%d]*%f*CWinv; ' % (
                species.id, species.id, species.weight) )

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckcpor(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get Cp/R as a function of T '))
        self._write(self.line('for all species (Eq 19)'))
        self._write('void CKCPOR'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpor)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('cp_R(cpor, tc);')
        
        self._outdent()

        self._write('}')

        return
    
    def _ckhort(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get H/RT as a function of T '))
        self._write(self.line('for all species (Eq 20)'))
        self._write('void CKHORT'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hort)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('speciesEnthalpy(hort, tc);')
        
        self._outdent()

        self._write('}')

        return
 
    def _cksor(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get S/R as a function of T '))
        self._write(self.line('for all species (Eq 21)'))
        self._write('void CKSOR'+sym+'(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sor)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('speciesEntropy(sor, tc);')
        
        self._outdent()

        self._write('}')

        return


    def _ckqc(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write()
        self._write()
        self._write(self.line('Returns the rate of progress for each reaction'))
        self._write('void CKQC'+sym+'(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        # convert C to SI units
        self._write()
        self._write(self.line('convert to SI'))
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e6;')
        self._outdent()
        self._write('}')
        
        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('progressRate(qdot, C, *T);')

        # convert C to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')

        # convert qdot to chemkin units
        self._write()
        self._write('for (id = 0; id < %d; ++id) {' % nReactions)
        self._indent()
        self._write('qdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return

    
    def _ckkfkr(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write()
        self._write(self.line('Returns the progress rates of each reactions'))
        self._write(self.line('Given P, T, and mole fractions'))
        self._write('void CKKFKR'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % nSpecies + self.line('temporary storage'))
        
        self._write('double PORT = 1e6 * (*P)/(%g * (*T)); ' % (R*kelvin*mole/erg) +
                    self.line('1e6 * P/RT so c goes to SI units'))
        
        # now compute conversion
        self._write()
        self._write(self.line('Compute conversion, see Eq 10'))
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('c[id] = x[id]*PORT;')
        self._outdent()
        self._write('}')
        
        # call progressRateFR
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('progressRateFR(q_f, q_r, c, *T);')

        # convert qdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nReactions )
        self._indent()
        self._write('q_f[id] *= 1.0e-6;')
        self._write('q_r[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckqyp(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write()
        self._write(self.line('Returns the progress rates of each reactions'))
        self._write(self.line('Given P, T, and mass fractions'))
        self._write('void CKQYP'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % nSpecies + self.line('temporary storage'))
        self._write('double YOW = 0; ')
        self._write('double PWORT; ')
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))
 
        self._write(self.line('PW/RT (see Eq. 7)'))
        self._write('PWORT = (*P)/(YOW * %g * (*T)); ' % (R*kelvin*mole/erg) )
        
        self._write(self.line('multiply by 1e6 so c goes to SI'))
        self._write('PWORT *= 1e6; ')

        # now compute conversion
        self._write(self.line('Now compute conversion (and go to SI)'))
        for species in self.species:
            self._write('c[%d] = PWORT * y[%d]*imw[%d]; ' % (
                species.id, species.id, species.id) )

        # call progressRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('progressRate(qdot, c, *T);')

        # convert qdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nReactions )
        self._indent()
        self._write('qdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckqxp(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write()
        self._write(self.line('Returns the progress rates of each reactions'))
        self._write(self.line('Given P, T, and mole fractions'))
        self._write('void CKQXP'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % nSpecies + self.line('temporary storage'))
        
        self._write('double PORT = 1e6 * (*P)/(%g * (*T)); ' % (R*kelvin*mole/erg) +
                    self.line('1e6 * P/RT so c goes to SI units'))
        
        # now compute conversion
        self._write()
        self._write(self.line('Compute conversion, see Eq 10'))
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('c[id] = x[id]*PORT;')
        self._outdent()
        self._write('}')
        
        # call progressRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('progressRate(qdot, c, *T);')

        # convert qdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nReactions )
        self._indent()
        self._write('qdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckqyr(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write()
        self._write(self.line('Returns the progress rates of each reactions'))
        self._write(self.line('Given rho, T, and mass fractions'))
        self._write('void CKQYR'+sym+'(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % nSpecies + self.line('temporary storage'))

        # now compute conversion
        self._write(self.line('See Eq 8 with an extra 1e6 so c goes to SI'))
        for species in self.species:
            self._write('c[%d] = 1e6 * (*rho) * y[%d]*imw[%d]; ' % (
                species.id, species.id, species.id) )
            
        # call progressRate
        self._write()
        self._write(self.line('call progressRate'))
        self._write('progressRate(qdot, c, *T);')

        # convert qdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nReactions )
        self._indent()
        self._write('qdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckqxr(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write()
        self._write(self.line('Returns the progress rates of each reactions'))
        self._write(self.line('Given rho, T, and mole fractions'))
        self._write('void CKQXR'+sym+'(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % nSpecies + self.line('temporary storage'))
        
        self._write('double XW = 0; '+self.line('See Eq 4, 11 in CK Manual'))
        self._write('double ROW; ')
        
        # compute mean molecular weight first (eq 3)
        self._write(self.line('Compute mean molecular wt first'))
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        # now compute conversion
        self._write(self.line('Extra 1e6 factor to take c to SI'))
        self._write('ROW = 1e6*(*rho) / XW;')
        self._write()
        self._write(self.line('Compute conversion, see Eq 11'))
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('c[id] = x[id]*ROW;')
        self._outdent()
        self._write('}')
        
        # call progressRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('progressRate(qdot, c, *T);')

        # convert qdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nReactions )
        self._indent()
        self._write('qdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return

    
    def __ckeqcontent(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double gort[%d]; ' % nSpecies + self.line(' temporary storage'))

        # compute the gibbs free energy
        self._write()
        self._write(self.line('compute the Gibbs free energy'))
        self._write('gibbs(gort, tc);')

        # compute the equilibrium constants
        self._write()
        self._write(self.line('compute the equilibrium constants'))
        self._write('equilibriumConstants(eqcon, gort, tT);')

        for reaction in mechanism.reaction():

            self._write()
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

            somepow = 0
            for symbol, coefficient in reaction.reactants:
                somepow = somepow - coefficient

            for symbol, coefficient in reaction.products:
                somepow = somepow + coefficient

            if somepow == 0:
                self._write(self.line(
                    'eqcon[%d] *= %g; ' % (reaction.id-1, (1e-6)**somepow) ) )
                
            else:
                self._write( 'eqcon[%d] *= %g; ' % (reaction.id-1, (1e-6)**somepow) ) 


    def _ckeqc(self, mechanism):

        self._write()
        self._write()
        self._write(self.line('Returns the equil constants for each reaction'))
        self._write('void CKEQC'+sym+'(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict eqcon)')
        self._write('{')
        self._indent()

        self.__ckeqcontent(mechanism)
                
        self._outdent()

        self._write('}')

        return

    
    def _ckeqyp(self, mechanism):

        import pyre
        periodic = pyre.handbook.periodicTable()

        self._write()
        self._write()
        self._write(self.line('Returns the equil constants for each reaction'))
        self._write(self.line('Given P, T, and mass fractions'))
        self._write('void CKEQYP'+sym+'(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)')
        self._write('{')
        self._indent()

        self.__ckeqcontent(mechanism)
        
        self._outdent()

        self._write('}')

        return


    def _ckeqxp(self, mechanism):

        import pyre
        periodic = pyre.handbook.periodicTable()

        self._write()
        self._write()
        self._write(self.line('Returns the equil constants for each reaction'))
        self._write(self.line('Given P, T, and mole fractions'))
        self._write('void CKEQXP'+sym+'(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)')
        self._write('{')
        self._indent()

        self.__ckeqcontent(mechanism)
        
        self._outdent()

        self._write('}')

        return


    def _ckeqyr(self, mechanism):

        import pyre
        periodic = pyre.handbook.periodicTable()

        self._write()
        self._write()
        self._write(self.line('Returns the equil constants for each reaction'))
        self._write(self.line('Given rho, T, and mass fractions'))
        self._write('void CKEQYR'+sym+'(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)')
        self._write('{')
        self._indent()

        self.__ckeqcontent(mechanism)
        
        self._outdent()

        self._write('}')

        return


    def _ckeqxr(self, mechanism):

        import pyre
        periodic = pyre.handbook.periodicTable()

        self._write()
        self._write()
        self._write(self.line('Returns the equil constants for each reaction'))
        self._write(self.line('Given rho, T, and mole fractions'))
        self._write('void CKEQXR'+sym+'(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)')
        self._write('{')
        self._indent()

        self.__ckeqcontent(mechanism)
        
        self._outdent()

        self._write('}')

        return


# Fuego Extensions. All functions in this section has the fe prefix
# All fuctions in this section uses the standard fuego chemkin functions
    def _ck_eytt(self, mechanism):

        nSpecies = len(mechanism.species())
        lowT,highT,dummy = self._analyzeThermodynamics(mechanism)
        
        self._write()
        self._write()
        self._write(self.line(
            'get temperature given internal energy in mass units and mass fracs'))
        self._write('int feeytt'+fsym+'(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t)')
        self._write('{')
        self._indent()

        self._write('const int maxiter = 50;')
        self._write('const double tol  = 0.001;')
        self._write('double ein  = *e;')
        self._write('double tmin = %g; // max lower bound for thermo def' % lowT)
        self._write('double tmax = %g; // min upper bound for thermo def' % highT)
        self._write('double e1,emin,emax,cv,t1,dt;')
        self._write('int i; // loop counter')
        self._write('CKUBMS'+sym+'(&tmin, y, iwrk, rwrk, &emin);')
        self._write('CKUBMS'+sym+'(&tmax, y, iwrk, rwrk, &emax);')
        self._write('if (ein < emin) {')
        self._indent()
        self._write(self.line('Linear Extrapolation below tmin'))
        self._write('CKCVBS'+sym+'(&tmin, y, iwrk, rwrk, &cv);')
        self._write('*t = tmin - (emin-ein)/cv;')
        self._write('return 1;')
        self._outdent()
        self._write('}')
        
        self._write('if (ein > emax) {')
        self._indent()
        self._write(self.line('Linear Extrapolation above tmax'))
        self._write('CKCVBS'+sym+'(&tmax, y, iwrk, rwrk, &cv);')
        self._write('*t = tmax - (emax-ein)/cv;')
        self._write('return 1;')
        self._outdent()
        self._write('}')

        self._write('t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);')
        self._write('for (i = 0; i < maxiter; ++i) {')
        self._indent()
        self._write('CKUBMS'+sym+'(&t1,y,iwrk,rwrk,&e1);')
        self._write('CKCVBS'+sym+'(&t1,y,iwrk,rwrk,&cv);')
        self._write('dt = (ein - e1) / cv;')
        self._write('if (dt > 100) { dt = 100; }')
        self._write('else if (dt < -100) { dt = -100; }')
        self._write('else if (fabs(dt) < tol) break;')
        self._write('t1 += dt;')
        self._outdent()
        self._write('}')
        
        self._write('*t = t1;')
        self._write('return 0;')
        
        self._outdent()

        self._write('}')

        return

 
    def _ck_hytt(self, mechanism):

        nSpecies = len(mechanism.species())
        lowT,highT,dummy = self._analyzeThermodynamics(mechanism)
        
        self._write()
        self._write()
        self._write(self.line(
            'get temperature given enthalpy in mass units and mass fracs'))
        self._write('int fehytt'+fsym+'(double * restrict h, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t)')
        self._write('{')
        self._indent()

        self._write('const int maxiter = 50;')
        self._write('const double tol  = 0.001;')
        self._write('double hin  = *h;')
        self._write('double tmin = %g; // max lower bound for thermo def' % lowT)
        self._write('double tmax = %g; // min upper bound for thermo def' % highT)
        self._write('double h1,hmin,hmax,cp,t1,dt;')
        self._write('int i; // loop counter')
        self._write('CKHBMS'+sym+'(&tmin, y, iwrk, rwrk, &hmin);')
        self._write('CKHBMS'+sym+'(&tmax, y, iwrk, rwrk, &hmax);')
        self._write('if (hin < hmin) {')
        self._indent()
        self._write(self.line('Linear Extrapolation below tmin'))
        self._write('CKCPBS'+sym+'(&tmin, y, iwrk, rwrk, &cp);')
        self._write('*t = tmin - (hmin-hin)/cp;')
        self._write('return 1;')
        self._outdent()
        self._write('}')
        
        self._write('if (hin > hmax) {')
        self._indent()
        self._write(self.line('Linear Extrapolation above tmax'))
        self._write('CKCPBS'+sym+'(&tmax, y, iwrk, rwrk, &cp);')
        self._write('*t = tmax - (hmax-hin)/cp;')
        self._write('return 1;')
        self._outdent()
        self._write('}')

        self._write('t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);')
        self._write('for (i = 0; i < maxiter; ++i) {')
        self._indent()
        self._write('CKHBMS'+sym+'(&t1,y,iwrk,rwrk,&h1);')
        self._write('CKCPBS'+sym+'(&t1,y,iwrk,rwrk,&cp);')
        self._write('dt = (hin - h1) / cp;')
        self._write('if (dt > 100) { dt = 100; }')
        self._write('else if (dt < -100) { dt = -100; }')
        self._write('else if (fabs(dt) < tol) break;')
        self._write('t1 += dt;')
        self._outdent()
        self._write('}')
        
        self._write('*t = t1;')
        self._write('return 0;')
        
        self._outdent()

        self._write('}')

        return

 
    def _ck_phity(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert phi[species] (specific mole nums) to y[species] (mass fracs)'))
        self._write('void fephity'+fsym+'(double * restrict phi, int * iwrk, double * restrict rwrk, double * restrict y)')
        self._write('{')
        self._indent()

        self._write('double XW  = 0; ')
        self._write('int id; ' + self.line('loop counter'))
        
        # compute mean molecular weight first (eq 3)
        self._write(self.line('Compute mean molecular wt first'))
        for species in self.species:
            self._write('y[%d] = phi[%d]*%f;   XW += y[%d]; ' % (
                species.id, species.id, species.weight, species.id) +
                        self.line('%s' % species.symbol))
 
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('y[id] = y[id]/XW;')
        self._outdent()
        self._write('}')
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
 
    def _ck_ytphi(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert y[species] (mass fracs) to phi[species] (specific mole num)'))
        self._write('void feytphi'+fsym+'(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict phi)')
        self._write('{')
        self._indent()

        for species in self.species:
            self._write('phi[%d] = y[%d]/%15.8e; ' % (
                species.id, species.id, species.weight/1000.0) +
                        self.line('%s (wt in kg)' % species.symbol))
 
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return


    def _ck_ctyr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'reverse of ytcr, useful for rate computations'))
        self._write('void fectyr'+fsym+'(double * restrict c, double * restrict rho, int * iwrk, double * restrict rwrk, double * restrict y)')
        self._write('{')
        self._indent()

        # now compute conversion
        for species in self.species:
            self._write('y[%d] = c[%d] * %f / (*rho); ' % (
                species.id, species.id, species.weight) )
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
 
    def _ck_cvrhs(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'ddebdf compatible right hand side of CV burner'))
        self._write(self.line(
            'rwrk[0] and rwrk[1] should contain rho and ene respectively'))
        self._write(self.line(
            'working variable phi contains specific mole numbers'))
        self._write('void fecvrhs'+fsym+'(double * restrict time, double * restrict phi, double * restrict phidot, double * restrict rwrk, int * iwrk)')

	self._write('{')
	self._indent()
	# main body
        self._write('double rho,ene; ' + self.line('CV Parameters'))
        self._write('double y[%s], wdot[%s]; ' % (self.nSpecies, self.nSpecies) +
                    self.line('temporary storage'))
        self._write('int i; ' + self.line('Loop counter'))
        self._write('double temperature,pressure; ' + self.line('temporary var'))
        self._write('rho = rwrk[0];')
        self._write('ene = rwrk[1];')
        self._write('fephity'+fsym+'(phi, iwrk, rwrk, y);')
        self._write('feeytt'+fsym+'(&ene, y, iwrk, rwrk, &temperature);')
        self._write('CKPY'+sym+'(&rho, &temperature,  y, iwrk, rwrk, &pressure);')
        self._write('CKWYP'+sym+'(&pressure, &temperature,  y, iwrk, rwrk, wdot);')
        self._write('for (i=0; i<%s; ++i) phidot[i] = wdot[i] / (rho/1000.0); ' % self.nSpecies)
        self._write()
        self._write('return;')

	self._outdent()
	self._write('}')
	return


    def _ck_cvdim(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'returns the dimensionality of the cv burner (number of species)'))
        self._write('int fecvdim'+fsym+'()')

	self._write('{')
	self._indent()
	# main body
        self._write('return %d;' % self.nSpecies)

	self._outdent()
	self._write('}')
	return

 
    def _ck_zndrhs(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'ddebdf compatible right hand side of ZND solver'))
        self._write(self.line( 'rwrk[0] : scaling factor for pressure'))
        self._write(self.line( 'rwrk[1] : preshock density (g/cc) '))
        self._write(self.line( 'rwrk[2] : detonation velocity (cm/s) '))
        self._write(self.line( 'solution vector: [P; rho; y0 ... ylast] '))
        self._write('void fezndrhs'+fsym+'(double * restrict time, double * restrict z, double * restrict zdot, double * restrict rwrk, int * iwrk)')

	self._write('{')
	self._indent()
	# main body
        self._write('double psc,rho1,udet; ' + self.line('ZND Parameters'))
        self._write('double wt[%s], hms[%s], wdot[%s]; ' %
                    (self.nSpecies, self.nSpecies, self.nSpecies) +
                    self.line('temporary storage'))
        self._write('int i; ' + self.line('Loop counter'))
        self._write(self.line('temporary variables'))
        self._write('double ru, T, uvel, wtm, p, rho, gam, son, xm, sum, drdy, eta, cp, cv ;')
        self._write('double * restrict y; ' + self.line('mass frac pointer'))
        self._write()
        self._write('ru = %g;' % (R * mole * kelvin / erg))
        self._write()
        self._write('psc = rwrk[0];')
        self._write('rho1 = rwrk[1];')
        self._write('udet = rwrk[2];')
        self._write()
        self._write('p = z[0] * psc;')
        self._write('rho = z[1];')
        self._write()
        self._write('y = &z[3];')
        self._write()
        self._write('CKMMWY'+sym+'(y, 0, 0, &wtm);')
        self._write()
        self._write('T = p * wtm / rho / ru;')
        self._write()
        self._write('uvel = (rho1 * udet)/ rho;')
        self._write()
        self._write('CKCPBS'+sym+'(&T, y, 0, 0, &cp);')
        self._write('CKCVBS'+sym+'(&T, y, 0, 0, &cv);')
        self._write('gam = cp/cv;')
        self._write()
        self._write('son = sqrt(fabs(gam*ru*T/wtm));')
        self._write('xm = uvel/son;')
        self._write()
        self._write('CKHMS'+sym+'(&T, 0, 0, hms);')
        self._write('CKWT'+sym+'(0, 0, wt);')
        self._write('CKWYP'+sym+'(&p, &T, y, 0, 0, wdot);')
        self._write()
        self._write('sum = 0.0;')
        self._write('for (i=0; i<%s; ++i) {' % self.nSpecies)
        self._indent()
        self._write('zdot[i+3] = wdot[i] * wt[i] / rho;')
        self._write('drdy = -rho * wtm / wt[i];')
        self._write('sum += -( drdy + rho * hms[i]/ (cp*T) ) * zdot[i+3];')
        self._outdent()
        self._write('}')
        self._write()
        self._write('eta = 1.0 - xm*xm;')
        self._write('zdot[0] = -(uvel*uvel/eta/psc)*sum;')
        self._write('zdot[1] = -sum/eta;')
        self._write('zdot[2] = uvel;')
        self._write()
        self._write('return;')

	self._outdent()
	self._write('}')
	return


    def _ck_znddim(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'returns the dimensionality of the ZND solver (3+number of species)'))
        self._write('int feznddim'+fsym+'()')

	self._write('{')
	self._indent()
	# main body
        self._write('return %d;' % (self.nSpecies + 3) )

	self._outdent()
	self._write('}')
	return
    
    def _ck_mechfile(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'returns the name of the source mechanism file '))
        self._write('char* femechfile'+fsym+'()')

	self._write('{')
	self._indent()
	# main body
        self._write('return "%s";' % mechanism.name())

	self._outdent()
	self._write('}')
	return

    def _ck_symnum(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'returns the species number'))
        self._write('int fesymnum'+fsym+'(const char* s1)')

	self._write('{')
	self._indent()
        
        for species in self.species:
            self._write('if (strcmp(s1, "%s")==0) return %d; ' % (
                species.symbol, species.id))
 
        self._write(self.line( 'species name not found' ))
        self._write('return -1;')

	self._outdent()
	self._write('}')
	return
    
    def _ck_symname(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'returns the species name'))
        self._write('char* fesymname'+fsym+'(int sn)')

	self._write('{')
	self._indent()

        for species in self.species:
            self._write('if (sn==%d) return "%s"; ' % (
                species.id, species.symbol))
 
        self._write(self.line( 'species name not found' ))
        self._write('return "NOTFOUND";')

	self._outdent()
	self._write('}')
	return
    
# Fuego's core routines section begins here
    def _molecularWeight(self, mechanism):

        import pyre
        periodic = pyre.handbook.periodicTable()
        
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        self._write(self.line('save molecular weights into array'))
        self._write('void molecularWeight(double * restrict wt)')
        self._write('{')
        self._indent()

        # molecular weights of all species
        for species in mechanism.species():

            weight = 0.0 #species.molecularWeight()
            for elem, coef in species.composition:
                aw = mechanism.element(elem).weight
                if not aw:
                    aw = periodic.symbol(elem.capitalize()).atomicWeight
                weight += coef * aw

            self._write('wt[%d] = %f; ' % (
                species.id, weight) + self.line('%s' % species.symbol))

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 

    def _atomicWeight(self, mechanism):

        self._write()
        self._write()
        self._write(self.line('save atomic weights into array'))
        self._write('void atomicWeight(double * restrict awt)')
        self._write('{')
        self._indent()
        import pyre
        periodic = pyre.handbook.periodicTable()
        for element in mechanism.element():
            aw = mechanism.element(element.symbol).weight
            if not aw:
                aw = periodic.symbol(element.symbol.capitalize()).atomicWeight

            self._write('awt[%d] = %f; ' % (
                element.id, aw) + self.line('%s' % element.symbol))

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 

    def _productionRate(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]

        if len(self.reactionIndex) != 7:
            print '\n\nCheck this!!!\n'
            sys.exit(1)
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        # OMP stuff
        self._write()
        self._write('static double T_save = -1;')
        self._write('#ifdef _OPENMP')
        self._write('#pragma omp threadprivate(T_save)')
        self._write('#endif')
        self._write()
        self._write('static double k_f_save[%d];' % nReactions)
        self._write('#ifdef _OPENMP')
        self._write('#pragma omp threadprivate(k_f_save)')
        self._write('#endif')
        self._write()
        self._write('static double Kc_save[%d];' % nReactions)
        self._write('#ifdef _OPENMP')
        self._write('#pragma omp threadprivate(Kc_save)')
        self._write('#endif')
        self._write()

        # main function
        self._write()
        self._write(self.line('compute the production rate for each species'))
        self._write('void productionRate(double * restrict wdot, double * restrict sc, double T)')
        self._write('{')
        self._indent()

        self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
        self._write('double invT = 1.0 / tc[1];')
        
        self._write()
        self._write('if (T != T_save)')
        self._write('{')
        self._indent()
        self._write('T_save = T;')
        self._write('comp_k_f(tc,invT,k_f_save);');
        self._write('comp_Kc(tc,invT,Kc_save);');
        self._outdent()
        self._write("}")

        self._write()
        self._write('double qdot, q_f[%d], q_r[%d];' % (nReactions,nReactions))
        self._write('comp_qfqr(q_f, q_r, sc, tc, invT);');

        self._write()
        self._write('for (int i = 0; i < %d; ++i) {' % nSpecies)
        self._indent()
        self._write('wdot[i] = 0.0;')
        self._outdent()
        self._write('}')

        for i in range(nReactions):
            self._write()
            self._write("qdot = q_f[%d]-q_r[%d];" % (i,i))
            reaction = mechanism.reaction(id=i)
            agents = list(set(reaction.reactants + reaction.products))
            agents = sorted(agents, key=lambda x: mechanism.species(x[0]).id)
            # note that a species might appear as both reactant and product
            # a species might alos appear twice or more on on each side
            # agents is a set that contains unique (symbol, coefficient)
            for a in agents:
                symbol, coefficient = a
                for b in reaction.reactants:
                    if b == a:
                        if coefficient == 1:
                            self._write("wdot[%d] -= qdot;" % (mechanism.species(symbol).id))
                        else:
                            self._write("wdot[%d] -= %d * qdot;" % (mechanism.species(symbol).id, coefficient))
                for b in reaction.products: 
                    if b == a:
                        if coefficient == 1:
                            self._write("wdot[%d] += qdot;" % (mechanism.species(symbol).id))
                        else:
                            self._write("wdot[%d] += %d * qdot;" % (mechanism.species(symbol).id, coefficient))


        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        # k_f function
        self._write()
        self._write('void comp_k_f(double * restrict tc, double invT, double * restrict k_f)')
        self._write('{')
        self._indent()
        self._outdent()
        self._write('#ifdef __INTEL_COMPILER')
        self._indent()
        self._write('#pragma simd')
        self._outdent()
        self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<%d; ++i) {' % (nReactions))
        self._indent()
        self._write("k_f[i] = prefactor_units[i] * fwd_A[i]")
        self._write("            * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);")
        self._outdent()
        self._write("};")
        self._write('return;')
        self._outdent()
        self._write('}')

        # Kc
        self._write()
        self._write('void comp_Kc(double * restrict tc, double invT, double * restrict Kc)')
        self._write('{')
        self._indent()

        self._write(self.line('compute the Gibbs free energy'))
        self._write('double g_RT[%d];' % (nSpecies))
        self._write('gibbs(g_RT, tc);')

        self._write()

        for reaction in mechanism.reaction():
            KcExpArg = self._sortedKcExpArg(mechanism, reaction)
            self._write("Kc[%d] = %s;" % (reaction.id-1,KcExpArg))

        self._write()
        
        self._outdent()
        self._write('#ifdef __INTEL_COMPILER')
        self._indent()
        self._write(' #pragma simd')
        self._outdent()
        self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<%d; ++i) {' % (nReactions))
        self._indent()
        self._write("Kc[i] = exp(Kc[i]);")
        self._outdent()
        self._write("};")

        self._write()

        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g * invT;' % (atm.value, R.value))
        self._write('double refCinv = 1 / refC;')

        self._write()

        for reaction in mechanism.reaction():
            KcConv = self._KcConv(mechanism, reaction)
            if KcConv:
                self._write("Kc[%d] *= %s;" % (reaction.id-1,KcConv))        
        
        self._write()

        self._write('return;')
        self._outdent()
        self._write('}')

        # qdot
        self._write()
        self._write('void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)')
        self._write('{')
        self._indent()

        nclassd = nReactions - nspecial
        nCorr   = n3body + ntroe + nsri + nlindemann

        for i in range(nclassd):
            self._write()
            reaction = mechanism.reaction(id=i)
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))
            self._write("qf[%d] = %s;" % (i, self._sortedPhaseSpace(mechanism, reaction.reactants)))
            if reaction.reversible:
                self._write("qr[%d] = %s;" % (i, self._sortedPhaseSpace(mechanism, reaction.products)))
            else:
                self._write("qr[%d] = 0.0;" % (i))
            if reaction.rev:
                print "reaction.rev not finished"
                sys.exit(1)

        self._write()
        self._write('double T = tc[1];')
        self._write()
        self._write(self.line('compute the mixture concentration'))
        self._write('double mixture = 0.0;')
        self._write('for (int i = 0; i < %d; ++i) {' % nSpecies)
        self._indent()
        self._write('mixture += sc[i];')
        self._outdent()
        self._write('}')

        self._write()
        self._write("double Corr[%d];" % nclassd)
        self._write('for (int i = 0; i < %d; ++i) {' % nclassd)
        self._indent()
        self._write('Corr[i] = 1.0;')
        self._outdent()
        self._write('}')

        if ntroe > 0:
            self._write()
            self._write(self.line(" troe"))
            self._write("{")
            self._indent()
            self._write("double alpha[%d];" % ntroe)
            alpha_d = {}
            for i in range(itroe[0],itroe[1]):
                ii = i - itroe[0]
                reaction = mechanism.reaction(id=i)
                if reaction.thirdBody:
                    alpha = self._enhancement(mechanism, reaction)
                    if alpha in alpha_d:
                        self._write("alpha[%d] = %s;" %(ii,alpha_d[alpha]))
                    else:
                        self._write("alpha[%d] = %s;" %(ii,alpha))
                        alpha_d[alpha] = "alpha[%d]" % ii

            if ntroe >= 4:
                self._outdent()
                self._outdent()
                self._write('#ifdef __INTEL_COMPILER')
                self._indent()
                self._indent()
                self._write(' #pragma simd')
                self._outdent()
                self._outdent()
                self._write('#endif')
                self._indent()
                self._indent()
            self._write("for (int i=%d; i<%d; i++)" %(itroe[0],itroe[1]))
            self._write("{")
            self._indent()
            self._write("double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;")
            self._write("redP = alpha[i-%d] / k_f_save[i] * phase_units[i] * low_A[i] * exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] *invT);" % itroe[0])
            self._write("F = redP / (1.0 + redP);")
            self._write("logPred = log10(redP);")
            self._write('logFcent = log10(')
            self._write('    (fabs(troe_Tsss[i]) > 1.e-100 ? (1.-troe_a[i])*exp(-T/troe_Tsss[i]) : 0.) ')
            self._write('    + (fabs(troe_Ts[i]) > 1.e-100 ? troe_a[i] * exp(-T/troe_Ts[i]) : 0.) ')
            self._write('    + (troe_len[i] == 4 ? exp(-troe_Tss[i] * invT) : 0.) );')
            self._write("troe_c = -.4 - .67 * logFcent;")
            self._write("troe_n = .75 - 1.27 * logFcent;")
            self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));")
            self._write("F_troe = pow(10., logFcent / (1.0 + troe*troe));")
            self._write("Corr[i] = F * F_troe;")
            self._outdent()
            self._write('}')

            self._outdent()
            self._write("}")

        if nsri > 0:
            self._write()
            self._write(self.line(" SRI"))
            self._write("{")
            self._indent()
            self._write("double alpha[%d];" % nsri)
            self._write("double redP, F, X, F_sri;")
            alpha_d = {}
            for i in range(isri[0],isri[1]):
                ii = i - isri[0]
                reaction = mechanism.reaction(id=i)
                if reaction.thirdBody:
                    alpha = self._enhancement(mechanism, reaction)
                    if alpha in alpha_d:
                        self._write("alpha[%d] = %s;" %(ii,alpha_d[alpha]))
                    else:
                        self._write("alpha[%d] = %s;" %(ii,alpha))
                        alpha_d[alpha] = "alpha[%d]" % ii

            if nsri >= 4:
                self._outdent()
                self._outdent()
                self._write('#ifdef __INTEL_COMPILER')
                self._indent()
                self._indent()
                self._write(' #pragma simd')
                self._outdent()
                self._outdent()
                self._write('#endif')
                self._indent()
                self._indent()
            self._write("for (int i=%d; i<%d; i++)" %(isri[0],isri[1]))
            self._write("{")
            self._indent()
            self._write("redP = alpha[i-%d] / k_f_save[i] * phase_units[i] * low_A[i] * exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] *invT);" % itroe[0])
            self._write("F = redP / (1.0 + redP);")
            self._write("logPred = log10(redP);")
            self._write("X = 1.0 / (1.0 + logPred*logPred);")
            self._write("F_sri = exp(X * log(sri_a[i] * exp(-sri_b[i]*invT)")
            self._write("   +  (sri_c[i] > 1.e-100 ? exp(T/sri_c[i]) : 0.0) )")
            self._write("   *  (sri_len[i] > 3 ? sri_d[i]*exp(sri_e[i]*tc[0]) : 1.0);")
            self._write("Corr[i] = F * F_sri;")
            self._outdent()
            self._write('}')

            self._outdent()
            self._write("}")

        if nlindemann > 0:
            self._write()
            self._write(self.line(" Lindemann"))
            self._write("{")
            self._indent()
            if nlindemann > 1:
                self._write("double alpha[%d];" % nlindemann)
            else:
                self._write("double alpha;")

            for i in range(ilindemann[0],ilindemann[1]):
                ii = i - ilindemann[0]
                reaction = mechanism.reaction(id=i)
                if reaction.thirdBody:
                    alpha = self._enhancement(mechanism, reaction)
                    if nlindemann > 1:
                        self._write("alpha[%d] = %s;" %(ii,alpha))
                    else:
                        self._write("alpha = %s;" %(alpha))

            if nlindemann == 1:
                self._write("double redP = alpha / k_f_save[%d] * phase_units[%d] * low_A[%d] * exp(low_beta[%d] * tc[0] - activation_units[%d] * low_Ea[%d] * invT);" 
                            % (ilindemann[0],ilindemann[0],ilindemann[0],ilindemann[0],ilindemann[0],ilindemann[0]))
                self._write("Corr[%d] = redP / (1. + redP);" % ilindemann[0])
            else:
                if nlindemann >= 4:
                    self._outdent()
                    self._write('#ifdef __INTEL_COMPILER')
                    self._indent()
                    self._write(' #pragma simd')
                    self._outdent()
                    self._write('#endif')
                    self._indent()
                self._write("for (int i=%d; i<%d; i++)" % (ilindemann[0], ilindemann[1]))
                self._write("{")
                self._indent()
                self._write("double redP = alpha[i-%d] / k_f_save[i] * phase_units[i] * low_A[i] * exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] * invT);"
                            % ilindemann[0])
                self._write("Corr[i] = redP / (1. + redP);")
                self._outdent()
                self._write('}')

            self._outdent()
            self._write("}")

        if n3body > 0:
            self._write()
            self._write(self.line(" simple three-body correction"))
            self._write("{")
            self._indent()
            self._write("double alpha;")
            alpha_save = ""
            for i in range(i3body[0],i3body[1]):
                reaction = mechanism.reaction(id=i)
                if reaction.thirdBody:
                    alpha = self._enhancement(mechanism, reaction)
                    if alpha != alpha_save:
                        alpha_save = alpha
                        self._write("alpha = %s;" % alpha)
                    self._write("Corr[%d] = alpha;" % i)
            self._outdent()
            self._write("}")

        self._write()
        self._write("for (int i=0; i<%d; i++)" % nclassd)
        self._write("{")
        self._indent()
        self._write("qf[i] *= Corr[i] * k_f_save[i];")
        self._write("qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];")
        self._outdent()
        self._write("}")
        
        if nspecial > 0:

            print "\n\n ***** WARNING: %d unclassified reactions\n" % nspecial

            self._write()
            self._write(self.line('unclassified reactions'))
            self._write('{')
            self._indent()

            self._write(self.line("reactions: %d to %d" % (ispecial[0]+1,ispecial[1])))

            self._write('double Kc;                      ' + self.line('equilibrium constant'))
            self._write('double k_f;                     ' + self.line('forward reaction rate'))
            self._write('double k_r;                     ' + self.line('reverse reaction rate'))
            self._write('double q_f;                     ' + self.line('forward progress rate'))
            self._write('double q_r;                     ' + self.line('reverse progress rate'))
            self._write('double phi_f;                   '
                        + self.line('forward phase space factor'))
            self._write('double phi_r;                   ' + self.line('reverse phase space factor'))
            self._write('double alpha;                   ' + self.line('enhancement'))

            self._write('double redP;                    ' + self.line('reduced pressure'))
            self._write('double logPred;                 ' + self.line('log of above'))
            self._write('double F;                       ' + self.line('fallof rate enhancement'))
            self._write()
            self._write('double F_troe;                  ' + self.line('TROE intermediate'))
            self._write('double logFcent;                ' + self.line('TROE intermediate'))
            self._write('double troe;                    ' + self.line('TROE intermediate'))
            self._write('double troe_c;                  ' + self.line('TROE intermediate'))
            self._write('double troe_n;                  ' + self.line('TROE intermediate'))

            for i in range(ispecial[0],ispecial[1]):
                self._write()
                reaction = mechanism.reaction(id=i)
                self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

                # compute the rates
                self._forwardRate(mechanism, reaction)
                self._reverseRate(mechanism, reaction)

                # store the progress rate
                self._write("qf[%d] = q_f;" % i)
                self._write("qr[%d] = q_r;" % i)

            self._outdent()
            self._write('}')

        self._write()
        self._write('return;')
        self._outdent()
        self._write('}')

        return


    def _DproductionRate(self, mechanism):

        species_list = [x.symbol for x in mechanism.species()]
        nSpecies = len(species_list)

        self._write()
        self._write(self.line('compute the reaction Jacobian'))
        self._write('void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.e6 * sc[k];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian(J, c, *Tp, *consP);')

        self._write()
        self._write('/* dwdot[k]/dT */')
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('J[%d+k] *= 1.e-6;' % (nSpecies*(nSpecies+1)))
        self._outdent()
        self._write('}')

        self._write()
        self._write('/* dTdot/d[X] */')
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('J[k*%d+%d] *= 1.e6;' % (nSpecies+1, nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return


    def _ajac(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write()
        self._write(self.line('compute the reaction Jacobian'))
        self._write('void aJacobian(double * restrict J, double * restrict sc, double T, int consP)')
        self._write('{')
        self._indent()

        self._write('for (int i=0; i<%d; i++) {' % (nSpecies+1)**2)
        self._indent()
        self._write('J[i] = 0.0;')
        self._outdent()
        self._write('}')
        
        self._write()

        self._write('double wdot[%d];' % (nSpecies))
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies))
        self._indent()
        self._write('wdot[k] = 0.0;')
        self._outdent()
        self._write('}')
        
        self._write()

        self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
        self._write('double invT = 1.0 / tc[1];')
        self._write('double invT2 = invT * invT;')

        self._write()

        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g / T;' % (atm.value, R.value))
        self._write('double refCinv = 1.0 / refC;')

        self._write()

        self._write(self.line('compute the mixture concentration'))
        self._write('double mixture = 0.0;')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('mixture += sc[k];')
        self._outdent()
        self._write('}')

        self._write()

        self._write(self.line('compute the Gibbs free energy'))
        self._write('double g_RT[%d];' % (nSpecies))
        self._write('gibbs(g_RT, tc);')

        self._write()

        self._write(self.line('compute the species enthalpy'))
        self._write('double h_RT[%d];' % (nSpecies))
        self._write('speciesEnthalpy(h_RT, tc);')

        self._write()

        self._write('double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;') 
        self._write('double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;')
        self._write('double dqdci, dcdc_fac, dqdc[%d];' % (nSpecies))
        self._write('double Pr, fPr, F, k_0, logPr;') 
        self._write('double logFcent, troe_c, troe_n, troePr_den, troePr, troe;')
        self._write('double Fcent1, Fcent2, Fcent3, Fcent;')
        self._write('double dlogFdc, dlogFdn, dlogFdcn_fac;')
        self._write('double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;')
        self._write('const double ln10 = log(10.0);')
        self._write('const double log10e = 1.0/log(10.0);')

        for i, reaction in zip(range(nReactions), mechanism.reaction()):

            lt = reaction.lt
            if lt:
                print "Landau-Teller reactions are not supported"
                sys.exit(1)

            self._write(self.line('reaction %d: %s' % (i+1, reaction.equation())))
            if reaction.low:  # case 1
                self._write(self.line('a pressure-fall-off reaction'))
                self._ajac_reaction(mechanism, reaction, 1)
            elif reaction.thirdBody:  # case 2
                self._write(self.line('a third-body and non-pressure-fall-off reaction'))
                self._ajac_reaction(mechanism, reaction, 2)
            else:  # case 3
                self._write(self.line('a non-third-body and non-pressure-fall-off reaction'))
                self._ajac_reaction(mechanism, reaction, 3)
            self._write()

        self._write('double c_R[%d], dcRdT[%d], e_RT[%d];' % (nSpecies, nSpecies, nSpecies))
        self._write('double * eh_RT;')
        self._write('if (consP) {')
        self._indent()

        self._write('cp_R(c_R, tc);')
        self._write('dcvpRdT(dcRdT, tc);')
        self._write('eh_RT = &h_RT[0];');

        self._outdent()
        self._write('}')
        self._write('else {')
        self._indent()

        self._write('cv_R(c_R, tc);')
        self._write('dcvpRdT(dcRdT, tc);')
        self._write('speciesInternalEnergy(e_RT, tc);');
        self._write('eh_RT = &e_RT[0];');

        self._outdent()
        self._write('}')

        self._write()

        self._write('double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('cmix += c_R[k]*sc[k];')
        self._write('dcmixdT += dcRdT[k]*sc[k];')
        self._write('ehmix += eh_RT[k]*wdot[k];')
        self._write('dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[%d+k];' % \
                        (nSpecies*(nSpecies+1)))
        self._outdent()
        self._write('}')

        self._write()
        self._write('double cmixinv = 1.0/cmix;')
        self._write('double tmp1 = ehmix*cmixinv;')
        self._write('double tmp3 = cmixinv*T;')
        self._write('double tmp2 = tmp1*tmp3;')
        self._write('double dehmixdc;')

        self._write('/* dTdot/d[X] */')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('dehmixdc = 0.0;')
        self._write('for (int m = 0; m < %d; ++m) {' % nSpecies)
        self._indent()
        self._write('dehmixdc += eh_RT[m]*J[k*%s+m];' % (nSpecies+1))
        self._outdent()
        self._write('}')        
        self._write('J[k*%d+%d] = tmp2*c_R[k] - tmp3*dehmixdc;' % (nSpecies+1,nSpecies))
        self._outdent()
        self._write('}')

        self._write('/* dTdot/dT */')
        self._write('J[%d] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;' % \
                        (nSpecies*(nSpecies+1)+nSpecies))

        self._outdent()
        self._write('}')
        return


    def _ajac_reaction(self, mechanism, reaction, rcase):

        if rcase == 1: # pressure-dependent reaction
            isPD = True
            if reaction.thirdBody:
                has_alpha = True
                self._write('/* also 3-body */')
            else:
                has_alpha = False
                self._write('/* non 3-body */')
                print 'FIXME: pressure dependent non-3-body reaction in _ajac_reaction'
                sys.exit(1)
        elif rcase == 2: # third-body and non-pressure-dependent reaction
            isPD = False
            has_alpha = True
        elif rcase == 3: # simple non-third and non-pressure-dependent reaction
            isPD = False
            has_alpha = False
        else:
            print '_ajac_reaction: wrong case ', rcase
            exit(1)

        nSpecies = len(mechanism.species())
        rea_dict = {}
        pro_dict = {}
        all_dict = {}
        sumNuk = 0
        for symbol, coefficient in reaction.reactants:
            k = mechanism.species(symbol).id
            sumNuk -= coefficient
            if k in rea_dict:
                coe_old = rea_dict[k][1]
                rea_dict[k] = (symbol, coefficient+coe_old)
            else:
                rea_dict[k] = (symbol,  coefficient)
        for symbol, coefficient in reaction.products:
            k = mechanism.species(symbol).id
            sumNuk += coefficient
            if k in pro_dict:
                coe_old = pro_dict[k][1]
                pro_dict[k] = (symbol, coefficient+coe_old)
            else:
                pro_dict[k] = (symbol, coefficient)
        for k in range(nSpecies):
            if k in rea_dict and k in pro_dict:
                sr, nur = rea_dict[k]
                sp, nup = pro_dict[k]
                all_dict[k] = (sr, nup-nur)
            elif k in rea_dict:
                sr, nur = rea_dict[k]
                all_dict[k] = (sr, -nur)
            elif k in pro_dict:
                sp, nup = pro_dict[k]
                all_dict[k] = (sp, nup)

        sorted_reactants = sorted(rea_dict.values())
        sorted_products = sorted(pro_dict.values()) 

        if not reaction.reversible:
            if isPD or has_alpha:
                print 'FIXME: inreversible reaction in _ajac_reaction may not work'
                self._write('/* FIXME: inreversible reaction in _ajac_reaction may not work*/')
            for k in range(nSpecies):
                if k in sorted_reactants and k in sorted_products:
                    print 'FIXME: inreversible reaction in _ajac_reaction may not work'
                    self._write('/* FIXME: inreversible reaction in _ajac_reaction may not work*/')

        if isPD:
            Corr_s = 'Corr *'
        elif has_alpha:
            Corr_s = 'alpha * '
        else:
            Corr_s = ''

        if has_alpha:
            self._write("/* 3-body correction factor */")
            self._write("alpha = %s;" % self._enhancement(mechanism, reaction))

        # forward
        self._write('/* forward */')
        self._write("phi_f = %s;" % self._sortedPhaseSpace(mechanism, sorted_reactants))
        #
        self._write("k_f = prefactor_units[%d] * fwd_A[%d]" % (reaction.id-1,reaction.id-1))
        self._write("            * exp(fwd_beta[%d] * tc[0] - activation_units[%d] * fwd_Ea[%d] * invT);"
                    %(reaction.id-1,reaction.id-1,reaction.id-1))
        self._write("dlnkfdT = fwd_beta[%d] * invT + activation_units[%d] * fwd_Ea[%d] * invT2;"
                    %(reaction.id-1,reaction.id-1,reaction.id-1))

        if isPD:
            self._write('/* pressure-fall-off */')
            self._write("k_0 = low_A[%d] * exp(low_beta[%d] * tc[0] - activation_units[%d] * low_Ea[%d] * invT);"
                        %(reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('Pr = phase_units[%d] * alpha / k_f * k_0;' % (reaction.id-1))
            self._write('fPr = Pr / (1.0+Pr);')
            self._write("dlnk0dT = low_beta[%d] * invT + activation_units[%d] * low_Ea[%d] * invT2;"
                        %(reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('dlogPrdT = log10e*(dlnk0dT - dlnkfdT);')
            self._write('dlogfPrdT = dlogPrdT / (1.0+Pr);')
            if reaction.sri:
                self._write('/* SRI form */')
                print "FIXME: sri not supported in _ajac_reaction yet"
                sys.exit(1)
            elif reaction.troe:
                self._write('/* Troe form */')
                troe = reaction.troe
                self._write("logPr = log10(Pr);")
                self._write('Fcent1 = (fabs(troe_Tsss[%d]) > 1.e-100 ? (1.-troe_a[%d])*exp(-T/troe_Tsss[%d]) : 0.);'
                            %(reaction.id-1,reaction.id-1,reaction.id-1))
                self._write('Fcent2 = (fabs(troe_Ts[%d]) > 1.e-100 ? troe_a[%d] * exp(-T/troe_Ts[%d]) : 0.);'
                            %(reaction.id-1,reaction.id-1,reaction.id-1))
                self._write('Fcent3 = (troe_len[%d] == 4 ? exp(-troe_Tss[%d] * invT) : 0.);'
                            %(reaction.id-1,reaction.id-1))
                self._write('Fcent = Fcent1 + Fcent2 + Fcent3;')
                self._write("logFcent = log10(Fcent);")
                self._write("troe_c = -.4 - .67 * logFcent;")
                self._write("troe_n = .75 - 1.27 * logFcent;")
                self._write("troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));")
                self._write("troePr = (troe_c + logPr) * troePr_den;")
                self._write("troe = 1.0 / (1.0 + troePr*troePr);")
                self._write("F = pow(10.0, logFcent * troe);")

                self._write("dlogFcentdT = log10e/Fcent*( ")
                self._write("    (fabs(troe_Tsss[%d]) > 1.e-100 ? -Fcent1/troe_Tsss[%d] : 0.)"
                            %(reaction.id-1,reaction.id-1))
                self._write("  + (fabs(troe_Ts[%d]) > 1.e-100 ? -Fcent2/troe_Ts[%d] : 0.)"
                            %(reaction.id-1,reaction.id-1))
                self._write("  + (troe_len[%d] == 4 ? Fcent3*troe_Tss[%d]*invT2 : 0.) );"
                            %(reaction.id-1,reaction.id-1))

                self._write("dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;")
                self._write('dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;')
                self._write('dlogFdn = dlogFdcn_fac * troePr;')
                self._write('dlogFdlogPr = dlogFdc;')
                self._write('dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;')
            else:
                self._write('/* Lindemann form */')
                self._write('F = 1.0;')
                self._write('dlogFdlogPr = 0.0;')
                self._write('dlogFdT = 0.0;')

        # reverse
        if not reaction.reversible:
            self._write('/* rate of progress */')
            if (not has_alpha) and (not isPD):
                self._write('q = k_f*phi_f;')
            else:
                self._write('q_nocor = k_f*phi_f;')
                if isPD:
                    self._write('Corr = fPr * F;')
                    self._write('q = Corr * q_nocor;')
                else: 
                    self._write('q = alpha * q_nocor;')

            if isPD:
                self._write('dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);')
                self._write('dqdT = %sdlnkfdT*k_f*phi_f + dlnCorrdT*q;' % Corr_s)
            else:
                self._write('dqdT = %sdlnkfdT*k_f*phi_f;' % Corr_s)
        else:
            self._write('/* reverse */')
            self._write("phi_r = %s;" % self._sortedPhaseSpace(mechanism, sorted_products))
            self._write('Kc = %s;' % self._sortedKc(mechanism, reaction))
            self._write('k_r = k_f / Kc;')

            dlnKcdT_s = 'invT * ('
            terms = []
            for symbol, coefficient in sorted(sorted_reactants,
                                              key=lambda x: mechanism.species(x[0]).id):
                k = mechanism.species(symbol).id
                if coefficient == 1:
                    terms.append('h_RT[%d]' % (k))
                else:
                    terms.append('%d*h_RT[%d]' % (coefficient, k))
            dlnKcdT_s += '-(' + ' + '.join(terms) + ')'
            terms = []
            for symbol, coefficient in sorted(sorted_products,
                                              key=lambda x: mechanism.species(x[0]).id):
                k = mechanism.species(symbol).id
                if coefficient == 1:
                    terms.append('h_RT[%d]' % (k))
                else:
                    terms.append('%d*h_RT[%d]' % (coefficient, k))
            dlnKcdT_s += ' + (' + ' + '.join(terms) + ')'
            if sumNuk > 0:
                dlnKcdT_s += ' - %d' % sumNuk
            elif sumNuk < 0:
                dlnKcdT_s += ' + %d' % (-sumNuk)
            dlnKcdT_s += ')'
            self._write('dlnKcdT = %s;' % dlnKcdT_s)

            self._write('dkrdT = (dlnkfdT - dlnKcdT)*k_r;')

            self._write('/* rate of progress */')
            if (not has_alpha) and (not isPD):
                self._write('q = k_f*phi_f - k_r*phi_r;')
            else:
                self._write('q_nocor = k_f*phi_f - k_r*phi_r;')
                if isPD:
                    self._write('Corr = fPr * F;')
                    self._write('q = Corr * q_nocor;')
                else: 
                    self._write('q = alpha * q_nocor;')

            if isPD:
                self._write('dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);')
                self._write('dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;' % Corr_s)
            else:
                self._write('dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r);' % Corr_s)

        self._write("/* update wdot */")
        for k in sorted(all_dict.keys()):
            s, nu = all_dict[k]
            if nu == 1:
                self._write('wdot[%d] += q; /* %s */' % (k, s))
            elif nu == -1:
                self._write('wdot[%d] -= q; /* %s */' % (k, s))
            elif nu > 0:
                self._write('wdot[%d] += %.17g * q; /* %s */' % (k, nu, s))
            elif nu < 0:
                self._write('wdot[%d] -= %.17g * q; /* %s */' % (k, -nu, s))

        if isPD:
            self._write('/* for convenience */')
            self._write('k_f *= Corr;')
            if reaction.reversible:
                self._write('k_r *= Corr;')
        elif has_alpha:
            self._write('/* for convenience */')
            self._write('k_f *= alpha;')
            if reaction.reversible:
                self._write('k_r *= alpha;')
            else:
                self._write('k_r = 0.0;')

        if isPD:
            self._write('dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);')

        def dqdc_simple(dqdc_s, k):
            if dqdc_s ==  "0":
                dqdc_s = ''
            if k in sorted(rea_dict.keys()):
                dps = self._DphaseSpace(mechanism,sorted_reactants,rea_dict[k][0])
                if dps == "1.0":
                    dps_s = ''
                else:
                    dps_s = '*'+dps
                dqdc_s += ' + k_f%s' % dps_s
            if reaction.reversible:
                if k in sorted(pro_dict.keys()):
                    dps = self._DphaseSpace(mechanism,sorted_products,pro_dict[k][0])
                    if dps == "1.0":
                        dps_s = ''
                    else:
                        dps_s = '*'+dps
                    dqdc_s += ' - k_r%s' % dps_s
            return dqdc_s

        if has_alpha or isPD:

            self._write('if (consP) {')
            self._indent()

            for k in range(nSpecies):
                dqdc_s = self._Denhancement(mechanism,reaction,k,True)
                if dqdc_s != "0":
                    if isPD:
                        if dqdc_s == "1":
                            dqdc_s ='dcdc_fac'
                        else:
                            dqdc_s +='*dcdc_fac'
                    elif has_alpha:
                        if dqdc_s == "1":
                            dqdc_s ='q_nocor'
                        else:
                            dqdc_s +='*q_nocor'

                dqdc_s = dqdc_simple(dqdc_s,k)
                if dqdc_s:
                    symb_k = self.species[k].symbol
                    self._write('/* d()/d[%s] */' % symb_k)
                    self._write('dqdci = %s;' % (dqdc_s))
                    #
                    for m in sorted(all_dict.keys()):
                        if all_dict[m][1] != 0:
                            s1 = 'J[%d] += %.17g * dqdci;' % (k*(nSpecies+1)+m, all_dict[m][1])
                            s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                            s2 = '/* dwdot[%s]/d[%s] */' % (all_dict[m][0], symb_k)
                            self._write(s1.ljust(30) + s2)

            self._outdent()
            self._write('}')
            self._write('else {')
            self._indent()

            for k in range(nSpecies):
                dqdc_s = self._Denhancement(mechanism,reaction,k,False)
                if dqdc_s != '0':
                    if dqdc_s == '1':
                        dqdc_s ='dcdc_fac'
                    elif isPD:
                        dqdc_s +='*dcdc_fac'
                elif has_alpha:
                    if dqdc_s == '1':
                        dqdc_s ='q_nocor'
                    else:
                        dqdc_s +='*q_nocor'
                dqdc_s = dqdc_simple(dqdc_s,k)
                if dqdc_s:
                    self._write('dqdc[%d] = %s;' % (k,dqdc_s))

            self._write('for (int k=0; k<%d; k++) {' % nSpecies)
            self._indent()
            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d*k+%d] += %.17g * dqdc[k];' % ((nSpecies+1), m, all_dict[m][1])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                    self._write(s1)
            self._outdent()
            self._write('}')

            self._outdent()
            self._write('}')

            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d] += %.17g * dqdT; /* dwdot[%s]/dT */' % \
                        (nSpecies*(nSpecies+1)+m, all_dict[m][1], all_dict[m][0])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                    self._write(s1)

        else:

            for k in range(nSpecies):
                dqdc_s = dqdc_simple('',k)
                if dqdc_s:
                    self._write('/* d()/d[%s] */' % all_dict[k][0])
                    self._write('dqdci = %s;' % (dqdc_s))
                    if reaction.reversible or k in rea_dict:
                        for m in sorted(all_dict.keys()):
                            if all_dict[m][1] != 0:
                                s1 = 'J[%d] += %.17g * dqdci;' % (k*(nSpecies+1)+m, all_dict[m][1])
                                s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                                s2 = '/* dwdot[%s]/d[%s] */' % (all_dict[m][0], all_dict[k][0])
                                self._write(s1.ljust(30) + s2)
            self._write('/* d()/dT */')
            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d] += %.17g * dqdT;' % (nSpecies*(nSpecies+1)+m, all_dict[m][1])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=').replace('+= -1 *', '-=')
                    s2 = '/* dwdot[%s]/dT */' % (all_dict[m][0])
                    self._write(s1.ljust(30) + s2)

        return


    def _vproductionRate(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        self._write()
        self._write()
        self._write(self.line('compute the production rate for each species'))
        self._write('void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)')
        self._write('{')
        self._indent()

        self._write('double k_f_s[%d*npt], Kc_s[%d*npt], mixture[npt], g_RT[%d*npt];'
                    % (nReactions, nReactions, nSpecies))
        self._write('double tc[5*npt], invT[npt];')

        self._write()

        self._outdent()
        self._write('#ifdef __INTEL_COMPILER')
        self._indent()
        self._write(' #pragma simd')
        self._outdent()
        self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        self._write('tc[0*npt+i] = log(T[i]);')
        self._write('tc[1*npt+i] = T[i];')
        self._write('tc[2*npt+i] = T[i]*T[i];')
        self._write('tc[3*npt+i] = T[i]*T[i]*T[i];')
        self._write('tc[4*npt+i] = T[i]*T[i]*T[i]*T[i];')
        self._write('invT[i] = 1.0 / T[i];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        self._write('mixture[i] = 0.0;')
        self._outdent()
        self._write('}')
        
        self._write()
        self._write('for (int n=0; n<%d; n++) {' % nSpecies)
        self._indent()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        self._write('mixture[i] += sc[n*npt+i];')
        self._write('wdot[n*npt+i] = 0.0;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()
        self._write('vcomp_k_f(npt, k_f_s, tc, invT);')
        self._write()
        self._write('vcomp_gibbs(npt, g_RT, tc);')
        self._write()
        self._write('vcomp_Kc(npt, Kc_s, g_RT, invT);')
        self._write()
        if nReactions <= 50:
            self._write('vcomp_wdot(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);')
        else:
            for i in range(0,nReactions,50):
                self._write('vcomp_wdot_%d_%d(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);' % (i+1,min(i+50,nReactions)))

        self._outdent()
        self._write('}')

        self._write()

        self._write('void vcomp_k_f(int npt, double * restrict k_f_s, double * restrict tc, double * restrict invT)')
        self._write('{')
        self._write('#ifdef __INTEL_COMPILER')
        self._indent()
        self._write('#pragma simd')
        self._outdent()
        self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        for reaction in mechanism.reaction():
            self._write("k_f_s[%d*npt+i] = prefactor_units[%d] * fwd_A[%d] * exp(fwd_beta[%d] * tc[i] - activation_units[%d] * fwd_Ea[%d] * invT[i]);" 
                        % (reaction.id-1,reaction.id-1,reaction.id-1,
                           reaction.id-1,reaction.id-1,reaction.id-1))
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')        

        self._write()

        self._write('void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)')
        self._write('{')
        self._indent()
        self._write(self.line('compute the Gibbs free energy'))
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        self._write('double tg[5], g[%d];' % nSpecies)
        self._write('tg[0] = tc[0*npt+i];')
        self._write('tg[1] = tc[1*npt+i];')
        self._write('tg[2] = tc[2*npt+i];')
        self._write('tg[3] = tc[3*npt+i];')
        self._write('tg[4] = tc[4*npt+i];')
        self._write()
        self._write('gibbs(g, tg);')
        self._write()
        for ispec in range(nSpecies):
            self._write('g_RT[%d*npt+i] = g[%d];' % (ispec, ispec))
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()

        self._write('void vcomp_Kc(int npt, double * restrict Kc_s, double * restrict g_RT, double * restrict invT)')
        self._write('{')
        self._write('#ifdef __INTEL_COMPILER')
        self._indent()
        self._write('#pragma simd')
        self._outdent()
        self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = (101325. / 8.31451) * invT[i];');
        self._write('double refCinv = 1.0 / refC;');
        self._write()
        for reaction in mechanism.reaction():
            K_c = self._vKc(mechanism, reaction)
            self._write("Kc_s[%d*npt+i] = %s;" % (reaction.id-1,K_c))
        self._outdent()
        self._write('}')        
        self._outdent()
        self._write('}')        

        self._write()
        if nReactions <= 50:
            self._write('void vcomp_wdot(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,')
            self._write('		double * restrict k_f_s, double * restrict Kc_s,')
            self._write('		double * restrict tc, double * restrict invT, double * restrict T)')
            self._write('{')
            self._vcomp_wdot(mechanism,0,nReactions)
            self._write('}')
        else:
            for i in range(0,nReactions,50):
                nr = min(50, nReactions-i)
                self._write('void vcomp_wdot_%d_%d(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,' % (i+1,i+nr))
                self._write('		double * restrict k_f_s, double * restrict Kc_s,')
                self._write('		double * restrict tc, double * restrict invT, double * restrict T)')
                self._write('{')
                self._vcomp_wdot(mechanism,i,nr)
                self._write('}')
                self._write()

        return

    def _vcomp_wdot(self, mechanism, istart, nr):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        self._write('#ifdef __INTEL_COMPILER')
        self._indent()
        self._write('#pragma simd')
        self._outdent()
        self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()

        self._write('double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;')
        if istart < isimple[0]:
            self._write('double alpha;')
        if istart < i3body[0]:
            self._write('double redP, F;') 
        if istart < ilindemann[0]:
            self._write('double logPred;')
            if ntroe>0:
                self._write('double logFcent, troe_c, troe_n, troe, F_troe;')
            if nsri>0:
                self._write('double X, F_sri;')

        first_id = istart + 1
        last_id  = istart + nr

        for reaction in mechanism.reaction():

            if reaction.id < first_id or reaction.id > last_id:
                continue

            self._write()
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

            # compute the rates
            self._vforwardRate(mechanism, reaction)
            self._vreverseRate(mechanism, reaction)

            # store the progress rate
            self._write("qdot = q_f - q_r;")

            agents = list(set(reaction.reactants + reaction.products))
            agents = sorted(agents, key=lambda x: mechanism.species(x[0]).id)
            # note that a species might appear as both reactant and product
            # a species might alos appear twice or more on on each side
            # agents is a set that contains unique (symbol, coefficient)
            for a in agents:
                symbol, coefficient = a
                for b in reaction.reactants:
                    if b == a:
                        if coefficient == 1:
                            self._write("wdot[%d*npt+i] -= qdot;" 
                                        % (mechanism.species(symbol).id))
                        else:
                            self._write("wdot[%d*npt+i] -= %d * qdot;" 
                                        % (mechanism.species(symbol).id, coefficient))
                for b in reaction.products: 
                    if b == a:
                        if coefficient == 1:
                            self._write("wdot[%d*npt+i] += qdot;" 
                                        % (mechanism.species(symbol).id))
                        else:
                            self._write("wdot[%d*npt+i] += %d * qdot;" 
                                        % (mechanism.species(symbol).id, coefficient))

        self._outdent()
        self._write('}')
        self._outdent()        

        return


    def _progressRate(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]

        if len(self.reactionIndex) != 7:
            print '\n\nCheck this!!!\n'
            sys.exit(1)
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        self._write()
        self._write()
        self._write(self.line('compute the progress rate for each reaction'))
        self._write('void progressRate(double * restrict qdot, double * restrict sc, double T)')
        self._write('{')
        self._indent()

        self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
        self._write('double invT = 1.0 / tc[1];')
        
        self._write()
        self._write('if (T != T_save)')
        self._write('{')
        self._indent()
        self._write('T_save = T;')
        self._write('comp_k_f(tc,invT,k_f_save);');
        self._write('comp_Kc(tc,invT,Kc_save);');
        self._outdent()
        self._write("}")

        self._write()
        self._write('double q_f[%d], q_r[%d];' % (nReactions,nReactions))
        self._write('comp_qfqr(q_f, q_r, sc, tc, invT);');

        self._write()
        self._write('for (int i = 0; i < %d; ++i) {' % nReactions)
        self._indent()
        self._write('qdot[i] = q_f[i] - q_r[i];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return


    def _initializeRateCalculation(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        # declarations
        self._write()
        self._write('int id; ' + self.line('loop counter'))

        self._write('double mixture;                 '
                    + self.line('mixture concentration'))
        self._write('double g_RT[%d];                ' % nSpecies
                    + self.line('Gibbs free energy'))

        self._write('double Kc;                      ' + self.line('equilibrium constant'))
        self._write('double k_f;                     ' + self.line('forward reaction rate'))
        self._write('double k_r;                     ' + self.line('reverse reaction rate'))
        self._write('double q_f;                     ' + self.line('forward progress rate'))
        self._write('double q_r;                     ' + self.line('reverse progress rate'))
        self._write('double phi_f;                   '
                    + self.line('forward phase space factor'))
        self._write('double phi_r;                   '
                    + self.line('reverse phase space factor'))
        self._write('double alpha;                   ' + self.line('enhancement'))


        self._write('double redP;                    ' + self.line('reduced pressure'))
        self._write('double logPred;                 ' + self.line('log of above'))
        self._write('double F;                       '
                    + self.line('fallof rate enhancement'))
        self._write()
        self._write('double F_troe;                  ' + self.line('TROE intermediate'))
        self._write('double logFcent;                ' + self.line('TROE intermediate'))
        self._write('double troe;                    ' + self.line('TROE intermediate'))
        self._write('double troe_c;                  ' + self.line('TROE intermediate'))
        self._write('double troe_n;                  ' + self.line('TROE intermediate'))
        self._write()

        self._write(
            'double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; '
            + self.line('temperature cache'))

        self._write()
        self._write('double invT = 1.0 / tc[1];')

        # compute the reference concentration
        self._write()
        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g / T;' % (atm.value, R.value))
        self._write('double refCinv = 1 / refC;')

        # compute the mixture concentration
        self._write()
        self._write(self.line('compute the mixture concentration'))
        self._write('mixture = 0.0;')
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('mixture += sc[id];')
        self._outdent()
        self._write('}')

        # compute the Gibbs free energies
        self._write()
        self._write(self.line('compute the Gibbs free energy'))
        self._write('gibbs(g_RT, tc);')
        
        return


    def _forwardRate(self, mechanism, reaction):

        lt = reaction.lt
        if lt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return self._landau(reaction)

        dim = self._phaseSpaceUnits(reaction.reactants)

        phi_f = self._phaseSpace(mechanism, reaction.reactants)
        self._write("phi_f = %s;" % phi_f)

        thirdBody = reaction.thirdBody
        if not thirdBody:
            self._write("k_f = k_f_save[%d];" % (reaction.id-1))
            self._write("q_f = phi_f * k_f;")
            return
            
        alpha = self._enhancement(mechanism, reaction)
        self._write("alpha = %s;" % alpha)

        sri = reaction.sri
        low = reaction.low
        troe = reaction.troe

        if not low:
            self._write("k_f = alpha * k_f_save[%d];" % (reaction.id-1))
            self._write("q_f = phi_f * k_f;")
            return

        self._write("k_f = k_f_save[%d];" % (reaction.id-1))

        self._write("redP = alpha / k_f * phase_units[%d] * low_A[%d] * exp(low_beta[%d] * tc[0] - activation_units[%d] * low_Ea[%d] *invT);"
                    %(reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))
        self._write("F = redP / (1 + redP);")

        if sri:
            self._write("logPred = log10(redP);")
            self._write("X = 1.0 / (1.0 + logPred*logPred);")
            self._write("F_sri = exp(X * log(sri_a[%d] * exp(-sri_b[%d]/T)" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   +  (sri_c[%d] > 1.e-100 ? exp(T/sri_c[%d]) : 0.) )" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   *  (sri_len[%d] > 3 ? sri_d[%d]*exp(sri_e[%d]*tc[0]) : 1);" 
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write("F *= F_sri;")

        elif troe:
            self._write("logPred = log10(redP);")

            self._write('logFcent = log10(')
            self._write('    (fabs(troe_Tsss[%d]) > 1.e-100 ? (1-troe_a[%d])*exp(-T/troe_Tsss[%d]) : 0) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (fabs(troe_Ts[%d]) > 1.e-100 ? troe_a[%d] * exp(-T/troe_Ts[%d]) : 0) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (troe_len[%d] == 4 ? exp(-troe_Tss[%d] * invT) : 0) );' 
                        % (reaction.id-1,reaction.id-1))
            self._write("troe_c = -.4 - .67 * logFcent;")
            self._write("troe_n = .75 - 1.27 * logFcent;")
            self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));")
            self._write("F_troe = pow(10, logFcent / (1.0 + troe*troe));")
            self._write("F *= F_troe;")

        self._write("k_f *= F;")
        self._write("q_f = phi_f * k_f;")
        return
        

    def _vforwardRate(self, mechanism, reaction):

        lt = reaction.lt
        if lt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return self._landau(reaction)

        dim = self._phaseSpaceUnits(reaction.reactants)

        phi_f = self._vphaseSpace(mechanism, reaction.reactants)
        self._write("phi_f = %s;" % phi_f)
                
        thirdBody = reaction.thirdBody
        if not thirdBody:
            self._write("k_f = k_f_s[%d*npt+i];" % (reaction.id-1))
            self._write("q_f = phi_f * k_f;")
            return
            
        alpha = self._venhancement(mechanism, reaction)
        self._write("alpha = %s;" % alpha)

        sri = reaction.sri
        low = reaction.low
        troe = reaction.troe

        if not low:
            self._write("k_f = alpha * k_f_s[%d*npt+i];" % (reaction.id-1))
            self._write("q_f = phi_f * k_f;")
            return

        self._write("k_f = k_f_s[%d*npt+i];" % (reaction.id-1))
        self._write("redP = alpha / k_f * phase_units[%d] * low_A[%d] * exp(low_beta[%d] * tc[i] - activation_units[%d] * low_Ea[%d] * invT[i]);" 
                    % (reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))
        self._write("F = redP / (1 + redP);")

        if sri:
            self._write("logPred = log10(redP);")
            self._write("X = 1.0 / (1.0 + logPred*logPred);")
            self._write("F_sri = exp(X * log(sri_a[%d] * exp(-sri_b[%d]/T[i])" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   +  (sri_c[%d] > 1.e-100 ? exp(T[i]/sri_c[%d]) : 0.) )" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   *  (sri_len[%d] > 3 ? sri_d[%d]*exp(sri_e[%d]*tc[i]) : 1.);" 
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write("F *= F_sri;")

        elif troe:
            self._write("logPred = log10(redP);")

            self._write('logFcent = log10(')
            self._write('    (fabs(troe_Tsss[%d]) > 1.e-100 ? (1.-troe_a[%d])*exp(-T[i]/troe_Tsss[%d]) : 0.) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (fabs(troe_Ts[%d]) > 1.e-100 ? troe_a[%d] * exp(-T[i]/troe_Ts[%d]) : 0.) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (troe_len[%d] == 4 ? exp(-troe_Tss[%d] * invT[i]) : 0.) );' 
                        % (reaction.id-1,reaction.id-1))
            
            d = .14
            self._write("troe_c = -.4 - .67 * logFcent;")
            self._write("troe_n = .75 - 1.27 * logFcent;")
            self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));")
            self._write("F_troe = pow(10., logFcent / (1.0 + troe*troe));")
            self._write("F *= F_troe;")

        self._write("k_f *= F;")
        self._write("q_f = phi_f * k_f;")
        return
        

    def _reverseRate(self, mechanism, reaction):
        if not reaction.reversible:
            self._write("q_r = 0.0;")
            return

        phi_r = self._phaseSpace(mechanism, reaction.products)
        self._write("phi_r = %s;" % phi_r)

        if reaction.rlt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return

        if reaction.rev:

            self._write("k_r = prefactor_units[%d] * rev_A[%d] * exp(rev_beta[%d] * tc[0] - activation_units[%d] * rev_Ea[%d] * invT);"
                        % (reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))

            thirdBody = reaction.thirdBody
            if thirdBody:
                self._write("k_r *= alpha;")

            self._write("q_r[%d] = phi_r * k_r;" % (reaction.id - 1))
            return
        
        self._write("Kc = Kc_save[%d];" % (reaction.id-1))

        self._write("k_r = k_f / Kc;")
        self._write("q_r = phi_r * k_r;")

        return


    def _vreverseRate(self, mechanism, reaction):
        if not reaction.reversible:
            self._write("q_r = 0.0;")
            return

        phi_r = self._vphaseSpace(mechanism, reaction.products)
        self._write("phi_r = %s;" % phi_r)

        if reaction.rlt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return

        if reaction.rev:

            self._write("k_r = prefactor_units[%d] * rev_A[%d] * exp(rev_beta[%d] * tc[i] - activation_units[%d] * rev_Ea[%d] * invT[i]);"
                        % (reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))

            thirdBody = reaction.thirdBody
            if thirdBody:
                self._write("k_r *= alpha;")

            self._write("q_f = phi_r * k_r;")
            return
        
        self._write("Kc = Kc_s[%d*npt+i];" % (reaction.id-1))

        self._write("k_r = k_f / Kc;")
        self._write("q_r = phi_r * k_r;")

        return


    def _progressRateFR(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]

        if len(self.reactionIndex) != 7:
            print '\n\nCheck this!!!\n'
            sys.exit(1)
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        self._write()
        self._write()
        self._write(self.line('compute the progress rate for each reaction'))
        self._write('void progressRateFR(double * restrict q_f, double * restrict q_r, double * restrict sc, double T)')
        self._write('{')
        self._indent()

        self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
        self._write('double invT = 1.0 / tc[1];')
        
        self._write()
        self._write('if (T != T_save)')
        self._write('{')
        self._indent()
        self._write('T_save = T;')
        self._write('comp_k_f(tc,invT,k_f_save);');
        self._write('comp_Kc(tc,invT,Kc_save);');
        self._outdent()
        self._write("}")

        self._write()
        self._write('comp_qfqr(q_f, q_r, sc, tc, invT);');

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return


    def _initializeRateCalculationFR(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        # declarations
        self._write()
        self._write('int id; ' + self.line('loop counter'))

        self._write('double mixture;                 '
                    + self.line('mixture concentration'))
        self._write('double g_RT[%d];                ' % nSpecies
                    + self.line('Gibbs free energy'))

        self._write('double Kc;                      ' + self.line('equilibrium constant'))
        self._write('double k_f;                     ' + self.line('forward reaction rate'))
        self._write('double k_r;                     ' + self.line('reverse reaction rate'))
        self._write('double phi_f;                   '
                    + self.line('forward phase space factor'))
        self._write('double phi_r;                   '
                    + self.line('reverse phase space factor'))
        self._write('double alpha;                   ' + self.line('enhancement'))


        self._write('double redP;                    ' + self.line('reduced pressure'))
        self._write('double logPred;                 ' + self.line('log of above'))
        self._write('double F;                       '
                    + self.line('fallof rate enhancement'))
        self._write()
        self._write('double F_troe;                  ' + self.line('TROE intermediate'))
        self._write('double logFcent;                ' + self.line('TROE intermediate'))
        self._write('double troe;                    ' + self.line('TROE intermediate'))
        self._write('double troe_c;                  ' + self.line('TROE intermediate'))
        self._write('double troe_n;                  ' + self.line('TROE intermediate'))
        self._write()

        self._write(
            'double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; '
            + self.line('temperature cache'))

        self._write()
        self._write('double invT = 1.0 / tc[1];')

        # compute the reference concentration
        self._write()
        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g / T;' % (atm.value, R.value))

        # compute the mixture concentration
        self._write()
        self._write(self.line('compute the mixture concentration'))
        self._write('mixture = 0.0;')
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('mixture += sc[id];')
        self._outdent()
        self._write('}')

        # compute the Gibbs free energies
        self._write()
        self._write(self.line('compute the Gibbs free energy'))
        self._write('gibbs(g_RT, tc);')
        
        return


    def _forwardRateFR(self, mechanism, reaction):

        lt = reaction.lt
        if lt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return self._landau(reaction)

        dim = self._phaseSpaceUnits(reaction.reactants)

        phi_f = self._phaseSpace(mechanism, reaction.reactants)
        self._write("phi_f = %s;" % phi_f)

        thirdBody = reaction.thirdBody
        if not thirdBody:
            self._write("k_f = k_f_save[%d];" % (reaction.id-1))
            self._write("q_f[%d] = phi_f * k_f;" % (reaction.id - 1))
            return
            
        alpha = self._enhancement(mechanism, reaction)
        self._write("alpha = %s;" % alpha)

        sri = reaction.sri
        low = reaction.low
        troe = reaction.troe

        if not low:
            self._write("k_f = alpha * k_f_save[%d];" % (reaction.id-1))
            self._write("q_f[%d] = phi_f * k_f;" % (reaction.id - 1))
            return

        self._write("k_f = k_f_save[%d];" % (reaction.id-1))

        self._write("redP = alpha / k_f * phase_units[%d] * low_A[%d] * exp(low_beta[%d] * tc[0] - activation_units[%d] * low_Ea[%d] *invT);"
                    % (reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))
        self._write("F = redP / (1 + redP);")

        if sri:
            self._write("logPred = log10(redP);")
            self._write("X = 1.0 / (1.0 + logPred*logPred);")
            self._write("F_sri = exp(X * log(sri_a[%d] * exp(-sri_b[%d]/T)" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   +  (sri_c[%d] > 1.e-100 ? exp(T/sri_c[%d]) : 0.) )" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   *  (sri_len[%d] > 3 ? sri_d[%d]*exp(sri_e[%d]*tc[0]) : 1.);" 
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write("F *= F_sri;")

        elif troe:
            self._write("logPred = log10(redP);")

            self._write('logFcent = log10(')
            self._write('    (fabs(troe_Tsss[%d]) > 1.e-100 ? (1.-troe_a[%d])*exp(-T/troe_Tsss[%d]) : 0.) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (fabs(troe_Ts[%d]) > 1.e-100 ? troe_a[%d] * exp(-T/troe_Ts[%d]) : 0.) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (troe_len[%d] == 4 ? exp(-troe_Tss[%d] * invT) : 0) );'
                        % (reaction.id-1,reaction.id-1))

            d = .14
            self._write("troe_c = -.4 - .67 * logFcent;")
            self._write("troe_n = .75 - 1.27 * logFcent;")
            self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));")
            self._write("F_troe = pow(10, logFcent / (1.0 + troe*troe));")
            self._write("F *= F_troe;")

        self._write("k_f *= F;")
        self._write("q_f[%d] = phi_f * k_f;" % (reaction.id - 1))
        return
        

    def _reverseRateFR(self, mechanism, reaction):
        if not reaction.reversible:
            self._write("q_r[%d] = 0.0;" % (reaction.id - 1))
            return

        phi_r = self._phaseSpace(mechanism, reaction.products)
        self._write("phi_r = %s;" % phi_r)

        if reaction.rlt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return

        if reaction.rev:
            self._write("k_r = prefactor_units[%d] * rev_A[%d] * exp(rev_beta[%d] * tc[0] - activation_units[%d] * rev_Ea[%d] * invT);"
                        % (reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))

            thirdBody = reaction.thirdBody
            if thirdBody:
                self._write("k_r *= alpha;")

            self._write("q_r[%d] = phi_r * k_r;" % (reaction.id - 1))
            return
        
        self._write("Kc = Kc_save[%d];" % (reaction.id-1))

        self._write("k_r = k_f / Kc;")

        self._write("q_r[%d] = phi_r * k_r;" % (reaction.id - 1))

        return

    def _prefactorUnits(self, code, exponent):

        if code == "mole/cm**3":
            units = mole / cm**3
        elif code == "moles":
            units = mole / cm**3
        elif code == "molecules":
            import pyre
            units = 1.0 / avogadro / cm**3
        else:
            import pyre
            pyre.debug.Firewall.hit("unknown prefactor units '%s'" % code)
            return 1

        return units ** exponent / second


    def _activationEnergyUnits(self, code):
        if code == "cal/mole":
            units = cal / mole
        elif code == "kcal/mole":
            units = kcal /mole
        elif code == "joules/mole":
            units = J / mole
        elif code == "kjoules/mole":
            units = kJ / mole
        elif code == "kelvins":
            units = Rc * kelvin
        else:
            pyre.debug.Firewall.hit("unknown activation energy units '%s'" % code)
            return 1

        return units


    def _equilibriumConstants(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('compute the equilibrium constants for each reaction'))
        self._write('void equilibriumConstants(double * restrict kc, double * restrict g_RT, double T)')
        self._write('{')
        self._indent()

        # compute the reference concentration
        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g / T;' % (atm.value, R.value))

        # compute the equilibrium constants
        for reaction in mechanism.reaction():
            self._write()
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

            K_c = self._Kc(mechanism, reaction)
            self._write("kc[%d] = %s;" % (reaction.id - 1, K_c))

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return


    def _phaseSpace(self, mechanism, reagents):

        phi = []

        for symbol, coefficient in reagents:
            conc = "sc[%d]" % mechanism.species(symbol).id
            phi += [conc] * coefficient

        return "*".join(phi)

    def _sortedPhaseSpace(self, mechanism, reagents):

        phi = []

        for symbol, coefficient in sorted(reagents,key=lambda x:mechanism.species(x[0]).id):
            conc = "sc[%d]" % mechanism.species(symbol).id
            phi += [conc] * coefficient

        return "*".join(phi)


    def _DphaseSpace(self, mechanism, reagents, r):

        phi = []

        for symbol, coefficient in sorted(reagents,key=lambda x:mechanism.species(x[0]).id):
            if symbol == r:
                if coefficient > 1:
                    conc = "sc[%d]" % mechanism.species(symbol).id
                    phi += ["%d" % coefficient]
                    phi += [conc] * (coefficient-1)
            else:
                conc = "sc[%d]" % mechanism.species(symbol).id
                phi += [conc] * coefficient

        if phi:
            return "*".join(phi)
        else:
            return "1.0"


    def _vphaseSpace(self, mechanism, reagents):

        phi = []

        for symbol, coefficient in sorted(reagents,key=lambda x:mechanism.species(x[0]).id):
            conc = "sc[%d*npt+i]" % mechanism.species(symbol).id
            phi += [conc] * coefficient

        return "*".join(phi)


    def _phaseSpaceUnits(self, reagents):
        dim = 0
        for symbol, coefficient in reagents:
            dim += coefficient

        return dim


    def _enhancement(self, mechanism, reaction):
        thirdBody = reaction.thirdBody
        if not thirdBody:
            import pyre
            pyre.debug.Firewall.hit("_enhancement called for a reaction without a third body")
            return

        species, coefficient = thirdBody
        efficiencies = reaction.efficiencies

        if not efficiencies:
            if species == "<mixture>":
                return "mixture"
            return "sc[%d]" % mechanism.species(species).id

        alpha = ["mixture"]
        for i, eff in enumerate(efficiencies):
            symbol, efficiency = eff
            factor = "(TB[%d][%d] - 1)" % (reaction.id-1, i)
            conc = "sc[%d]" % mechanism.species(symbol).id
            alpha.append("%s*%s" % (factor, conc))

        return " + ".join(alpha).replace('+ -','- ')

    def _Denhancement(self, mechanism, reaction, kid, consP):
        thirdBody = reaction.thirdBody
        if not thirdBody:
            import pyre
            pyre.debug.Firewall.hit("_enhancement called for a reaction without a third body")
            return

        species, coefficient = thirdBody
        efficiencies = reaction.efficiencies

        if not efficiencies:
            if species == "<mixture>":
                if consP:
                    return "0"
                else:
                    return "1"
            elif mechanism.species(species).id == kid:
                return "1"
            else:
                return "0"
        else:
            if consP:
                for i, eff in enumerate(efficiencies):
                    symbol, efficiency = eff
                    if mechanism.species(symbol).id == kid:
                        return "(TB[%d][%d] - 1)" % (reaction.id-1, i)
                return "0"
            else:
                for i, eff in enumerate(efficiencies):
                    symbol, efficiency = eff
                    if mechanism.species(symbol).id == kid:
                        return "TB[%d][%d]" % (reaction.id-1,i)
                return "1"

    def _venhancement(self, mechanism, reaction):
        thirdBody = reaction.thirdBody
        if not thirdBody:
            import pyre
            pyre.debug.Firewall.hit("_enhancement called for a reaction without a third body")
            return

        species, coefficient = thirdBody
        efficiencies = reaction.efficiencies

        if not efficiencies:
            if species == "<mixture>":
                return "mixture[i]"
            return "sc[%d*npt+i]" % mechanism.species(species).id

        alpha = ["mixture[i]"]
        for i, eff in enumerate(efficiencies):
            symbol, efficiency = eff
            factor = "(TB[%d][%d] - 1)" % (reaction.id-1,i)
            conc = "sc[%d*npt+i]" % mechanism.species(symbol).id
            alpha.append("%s*%s" % (factor, conc))

        return " + ".join(alpha)


    def _cv(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute Cv/R at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("cv_R", self._cvNASA, speciesInfo)

        return
    
    def _cp(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute Cp/R at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("cp_R", self._cpNASA, speciesInfo)

        return

    def _dcvpdT(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("dcvpRdT", self._dcpdTNASA, speciesInfo)

        return

    def _gibbs(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the g/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("gibbs", self._gibbsNASA, speciesInfo, 1)

        return

    def _helmholtz(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the a/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("helmholtz", self._helmholtzNASA, speciesInfo, 1)

        return

    def _speciesEntropy(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the S/R at the given temperature (Eq 21)'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("speciesEntropy", self._entropyNASA, speciesInfo)

        return

    def _speciesInternalEnergy(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the e/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("speciesInternalEnergy", self._internalEnergy, speciesInfo, 1)

        return

    def _speciesEnthalpy(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the h/(RT) at the given temperature (Eq 20)'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("speciesEnthalpy", self._enthalpyNASA, speciesInfo, 1)

        return

    

    def _generateThermoRoutine(self, name, expressionGenerator, speciesInfo, needsInvT=0):

        lowT, highT, midpoints = speciesInfo
        
        self._write('void %s(double * restrict species, double * restrict tc)' % name)
        self._write('{')

        self._indent()

        # declarations
        self._write()
        self._write(self.line('temperature'))
        self._write('double T = tc[1];')
        if needsInvT != 0:
           self._write('double invT = 1 / T;')
        if needsInvT == 2:
           self._write('double invT2 = invT*invT;')

        # temperature check
        # self._write()
        # self._write(self.line('check the temperature value'))
        # self._write('if (T < %g || T > %g) {' % (lowT, highT))
        # self._indent()
        # self._write(
        #     'fprintf(stderr, "temperature %%g is outside the range (%g, %g)", T);'
        #     % (lowT, highT))
        # self._write('return;')
        # self._outdent()
        # self._write('}')
                    
        for midT, speciesList in midpoints.items():

            self._write('')
            self._write(self.line('species with midpoint at T=%g kelvin' % midT))
            self._write('if (T < %g) {' % midT)
            self._indent()

            for species, lowRange, highRange in speciesList:
                self._write(self.line('species %d: %s' % (species.id, species.symbol)))
                self._write('species[%d] =' % species.id)
                self._indent()
                expressionGenerator(lowRange.parameters)
                self._outdent()

            self._outdent()
            self._write('} else {')
            self._indent()

            for species, lowRange, highRange in speciesList:
                self._write(self.line('species %d: %s' % (species.id, species.symbol)))
                self._write('species[%d] =' % species.id)
                self._indent()
                expressionGenerator(highRange.parameters)
                self._outdent()

            self._outdent()
            self._write('}')
            
        self._write('return;')
        self._outdent()

        self._write('}')

        return


    def _analyzeThermodynamics(self, mechanism):
        lowT = 0.0
        highT = 1000000.0

        midpoints = {}

        for species in mechanism.species():

            models = species.thermo
            if len(models) > 2:
                print 'species: ', species
                import pyre
                pyre.debug.Firewall.hit("unsupported configuration in species.thermo")
                return
            
            m1 = models[0]
            m2 = models[1]

            if m1.lowT < m2.lowT:
                lowRange = m1
                highRange = m2
            else:
                lowRange = m2
                highRange = m1

            low = lowRange.lowT
            mid = lowRange.highT
            high = highRange.highT

            if low > lowT:
                lowT = low
            if high < highT:
                highT = high

            midpoints.setdefault(mid, []).append((species, lowRange, highRange))
            
        return lowT, highT, midpoints


    def _Kc(self, mechanism, reaction):

        dim = 0
        dG = ""

        terms = []
        for symbol, coefficient in reaction.reactants:
            if coefficient == 1:
                factor = ""
            else:
                factor = "%d * " % coefficient
                    
            terms.append("%sg_RT[%d]" % (factor, mechanism.species(symbol).id))
            dim -= coefficient
        dG += '(' + ' + '.join(terms) + ')'

        # flip the signs
        terms = []
        for symbol, coefficient in reaction.products:
            if coefficient == 1:
                factor = ""
            else:
                factor = "%d * " % coefficient
            terms.append("%sg_RT[%d]" % (factor, mechanism.species(symbol).id))
            dim += coefficient
        dG += ' - (' + ' + '.join(terms) + ')'

        K_p = 'exp(' + dG + ')'

        if dim == 0:
            conversion = ""
        elif dim > 0:
            conversion = "*".join(["refC"] * dim) + ' * '
        else:
            conversion = "1.0 / (" + "*".join(["refC"] * abs(dim)) + ') * '

        K_c = conversion + K_p

        return K_c

    def _KcConv(self, mechanism, reaction):
        dim = 0
        for symbol, coefficient in reaction.reactants:
            dim -= coefficient
        # flip the signs
        for symbol, coefficient in reaction.products:
            dim += coefficient

        if dim == 0:
            conversion = ""
        elif dim > 0:
            conversion = "*".join(["refC"] * dim)
        else:
            conversion = "*".join(["refCinv"] * abs(dim))

        return conversion

    def _sortedKcExpArg(self, mechanism, reaction):

        nSpecies = len(mechanism.species())

        terms = []
        for i in range(nSpecies):
            terms.append('')
        for symbol, coefficient in reaction.reactants:
            if coefficient == 1:
                factor = " + "
            else:
                factor = " + %d*" % coefficient
            i = mechanism.species(symbol).id
            terms[i] += "%sg_RT[%d]"%(factor,i)

        for symbol, coefficient in reaction.products:
            if coefficient == 1:
                factor = " - "    # flip the signs
            else:
                factor = " - %d*" % coefficient
            i = mechanism.species(symbol).id
            terms[i] += "%sg_RT[%d]"%(factor,i)

        dG = ""
        for i in range(nSpecies):
            if terms[i]:
                dG += terms[i]
        if dG[0:3] == " + ":
            return dG[3:]
        else:
            return "-"+dG[3:]


    def _sortedKc(self, mechanism, reaction):
        conv = self._KcConv(mechanism,reaction)
        exparg = self._sortedKcExpArg(mechanism,reaction)
        if conv:
            return conv + ' * exp('+exparg+')'
        else:
            return 'exp('+exparg+')'


    def _vKc(self, mechanism, reaction):

        dim = 0
        dG = ""

        terms = []
        for symbol, coefficient in sorted(reaction.reactants, 
                                          key=lambda x: mechanism.species(x[0]).id):
            if coefficient == 1:
                factor = ""
            else:
                factor = "%d * " % coefficient
                    
            terms.append("%sg_RT[%d*npt+i]" % (factor, mechanism.species(symbol).id))
            dim -= coefficient
        dG += '(' + ' + '.join(terms) + ')'

        # flip the signs
        terms = []
        for symbol, coefficient in sorted(reaction.products,
                                          key=lambda x: mechanism.species(x[0]).id):
            if coefficient == 1:
                factor = ""
            else:
                factor = "%d * " % coefficient
            terms.append("%sg_RT[%d*npt+i]" % (factor, mechanism.species(symbol).id))
            dim += coefficient
        dG += ' - (' + ' + '.join(terms) + ')'

        K_p = 'exp(' + dG + ')'

        if dim == 0:
            conversion = ""
        elif dim > 0:
            conversion = "*".join(["refC"] * dim) + ' * '
        else:
            conversion = "*".join(["refCinv"] * abs(dim)) + ' * '

        K_c = conversion + K_p

        return K_c


    def _Kc_exparg(self, mechanism, reaction):

        dG = ""

        terms = []
        for symbol, coefficient in reaction.reactants:
            if coefficient == 1:
                factor = ""
            else:
                factor = "%d * " % coefficient
                    
            terms.append("%sg_RT[%d]" % (factor, mechanism.species(symbol).id))
        dG += '(' + ' + '.join(terms) + ')'

        # flip the signs
        terms = []
        for symbol, coefficient in reaction.products:
            if coefficient == 1:
                factor = ""
            else:
                factor = "%d * " % coefficient
            terms.append("%sg_RT[%d]" % (factor, mechanism.species(symbol).id))
        dG += ' - (' + ' + '.join(terms) + ')'

        K_p = 'exp(' + dG + ')'

        return dG

    def _cpNASA(self, parameters):
        self._write('%+15.8e' % parameters[0])
        self._write('%+15.8e * tc[1]' % parameters[1])
        self._write('%+15.8e * tc[2]' % parameters[2])
        self._write('%+15.8e * tc[3]' % parameters[3])
        self._write('%+15.8e * tc[4];' % parameters[4])
        return

    def _dcpdTNASA(self, parameters):
        self._write('%+15.8e' % parameters[1])
        self._write('%+15.8e * tc[1]' % (parameters[2]*2.))
        self._write('%+15.8e * tc[2]' % (parameters[3]*3.))
        self._write('%+15.8e * tc[3];' % (parameters[4]*4.))
        return

    def _cvNASA(self, parameters):
        self._write('%+15.8e' % (parameters[0] - 1.0))
        self._write('%+15.8e * tc[1]' % parameters[1])
        self._write('%+15.8e * tc[2]' % parameters[2])
        self._write('%+15.8e * tc[3]' % parameters[3])
        self._write('%+15.8e * tc[4];' % parameters[4])
        return


    def _enthalpyNASA(self, parameters):
        self._write('%+15.8e' % parameters[0])
        self._write('%+15.8e * tc[1]' % (parameters[1]/2))
        self._write('%+15.8e * tc[2]' % (parameters[2]/3))
        self._write('%+15.8e * tc[3]' % (parameters[3]/4))
        self._write('%+15.8e * tc[4]' % (parameters[4]/5))
        self._write('%+15.8e * invT;' % (parameters[5]))
        return


    def _internalEnergy(self, parameters):
        self._write('%+15.8e' % (parameters[0] - 1.0))
        self._write('%+15.8e * tc[1]' % (parameters[1]/2))
        self._write('%+15.8e * tc[2]' % (parameters[2]/3))
        self._write('%+15.8e * tc[3]' % (parameters[3]/4))
        self._write('%+15.8e * tc[4]' % (parameters[4]/5))
        self._write('%+15.8e * invT;' % (parameters[5]))
        return

    
    def _gibbsNASA(self, parameters):
        self._write('%+20.15e * invT' % parameters[5])
        self._write('%+20.15e' % (parameters[0] - parameters[6]))
        self._write('%+20.15e * tc[0]' % (-parameters[0]))
        self._write('%+20.15e * tc[1]' % (-parameters[1]/2))
        self._write('%+20.15e * tc[2]' % (-parameters[2]/6))
        self._write('%+20.15e * tc[3]' % (-parameters[3]/12))
        self._write('%+20.15e * tc[4];' % (-parameters[4]/20))
        return
    
    def _helmholtzNASA(self, parameters):
        self._write('%+15.8e * invT' % parameters[5])
        self._write('%+15.8e' % (parameters[0] - parameters[6] - 1.0))
        self._write('%+15.8e * tc[0]' % (-parameters[0]))
        self._write('%+15.8e * tc[1]' % (-parameters[1]/2))
        self._write('%+15.8e * tc[2]' % (-parameters[2]/6))
        self._write('%+15.8e * tc[3]' % (-parameters[3]/12))
        self._write('%+15.8e * tc[4];' % (-parameters[4]/20))
        return

    def _entropyNASA(self, parameters):
        self._write('%+15.8e * tc[0]' % parameters[0])
        self._write('%+15.8e * tc[1]' % (parameters[1]))
        self._write('%+15.8e * tc[2]' % (parameters[2]/2))
        self._write('%+15.8e * tc[3]' % (parameters[3]/3))
        self._write('%+15.8e * tc[4]' % (parameters[4]/4))
        self._write('%+15.8e ;' % (parameters[6]))
        return

    def _T_given_ey(self, mechanism):
        self._write(self.line(' get temperature given internal energy in mass units and mass fracs'))
        self._write('void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int * ierr)')
        self._write('{')
        self._write('#ifdef CONVERGENCE')
        self._indent()
        self._write('const int maxiter = 5000;')
        self._write('const double tol  = 1.e-12;')
        self._outdent()
        self._write('#else')
        self._indent()
        self._write('const int maxiter = 200;')
        self._write('const double tol  = 1.e-6;')
        self._outdent()
        self._write('#endif')
        self._indent()
        self._write('double ein  = *e;')
        self._write('double tmin = 250;'+self.line('max lower bound for thermo def'))
        self._write('double tmax = 4000;'+self.line('min upper bound for thermo def'))
        self._write('double e1,emin,emax,cv,t1,dt;')
        self._write('int i;'+self.line(' loop counter'))
        self._write('CKUBMS(&tmin, y, iwrk, rwrk, &emin);')
        self._write('CKUBMS(&tmax, y, iwrk, rwrk, &emax);')
        self._write('if (ein < emin) {')
        self._indent()
        self._write(self.line('Linear Extrapolation below tmin'))
        self._write('CKCVBS(&tmin, y, iwrk, rwrk, &cv);')
        self._write('*t = tmin - (emin-ein)/cv;')
        self._write('*ierr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('if (ein > emax) {')
        self._indent()
        self._write(self.line('Linear Extrapolation above tmax'))
        self._write('CKCVBS(&tmax, y, iwrk, rwrk, &cv);')
        self._write('*t = tmax - (emax-ein)/cv;')
        self._write('*ierr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('t1 = *t;')
        self._write('if (t1 < tmin || t1 > tmax) {')
        self._indent()
        self._write('t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);')
        self._outdent()
        self._write('}')
        self._write('for (i = 0; i < maxiter; ++i) {')
        self._indent()
        self._write('CKUBMS(&t1,y,iwrk,rwrk,&e1);')
        self._write('CKCVBS(&t1,y,iwrk,rwrk,&cv);')
        self._write('dt = (ein - e1) / cv;')
        self._write('if (dt > 100.) { dt = 100.; }')
        self._write('else if (dt < -100.) { dt = -100.; }')
        self._write('else if (fabs(dt) < tol) break;')
        self._write('else if (t1+dt == t1) break;')
        self._write('t1 += dt;')
        self._outdent()
        self._write('}')
        self._write('*t = t1;')
        self._write('*ierr = 0;')
        self._write('return;')
        self._outdent()
        self._write('}')

    def _T_given_hy(self, mechanism):
        self._write(self.line(' get temperature given enthalpy in mass units and mass fracs'))
        self._write('void GET_T_GIVEN_HY(double * restrict h, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int * ierr)')
        self._write('{')
        self._write('#ifdef CONVERGENCE')
        self._indent()
        self._write('const int maxiter = 5000;')
        self._write('const double tol  = 1.e-12;')
        self._outdent()
        self._write('#else')
        self._indent()
        self._write('const int maxiter = 200;')
        self._write('const double tol  = 1.e-6;')
        self._outdent()
        self._write('#endif')
        self._indent()
        self._write('double hin  = *h;')
        self._write('double tmin = 250;'+self.line('max lower bound for thermo def'))
        self._write('double tmax = 4000;'+self.line('min upper bound for thermo def'))
        self._write('double h1,hmin,hmax,cp,t1,dt;')
        self._write('int i;'+self.line(' loop counter'))
        self._write('CKHBMS(&tmin, y, iwrk, rwrk, &hmin);')
        self._write('CKHBMS(&tmax, y, iwrk, rwrk, &hmax);')
        self._write('if (hin < hmin) {')
        self._indent()
        self._write(self.line('Linear Extrapolation below tmin'))
        self._write('CKCPBS(&tmin, y, iwrk, rwrk, &cp);')
        self._write('*t = tmin - (hmin-hin)/cp;')
        self._write('*ierr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('if (hin > hmax) {')
        self._indent()
        self._write(self.line('Linear Extrapolation above tmax'))
        self._write('CKCPBS(&tmax, y, iwrk, rwrk, &cp);')
        self._write('*t = tmax - (hmax-hin)/cp;')
        self._write('*ierr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('t1 = *t;')
        self._write('if (t1 < tmin || t1 > tmax) {')
        self._indent()
        self._write('t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);')
        self._outdent()
        self._write('}')
        self._write('for (i = 0; i < maxiter; ++i) {')
        self._indent()
        self._write('CKHBMS(&t1,y,iwrk,rwrk,&h1);')
        self._write('CKCPBS(&t1,y,iwrk,rwrk,&cp);')
        self._write('dt = (hin - h1) / cp;')
        self._write('if (dt > 100.) { dt = 100.; }')
        self._write('else if (dt < -100.) { dt = -100.; }')
        self._write('else if (fabs(dt) < tol) break;')
        self._write('else if (t1+dt == t1) break;')
        self._write('t1 += dt;')
        self._outdent()
        self._write('}')
        self._write('*t = t1;')
        self._write('*ierr = 0;')
        self._write('return;')
        self._outdent()
        self._write('}')
# version
__id__ = "$Id$"

#  End of file 
