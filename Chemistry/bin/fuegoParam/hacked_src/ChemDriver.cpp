#include <winstd.H>
#include <iostream>

#include "ChemDriver.H"
#include "ChemDriver_F.H"
#include <ParallelDescriptor.H>
#include <ParmParse.H>

const Real HtoTerrMAX_DEF  = 1.e-8;
const int  HtoTiterMAX_DEF = 20;
const Real Tmin_trans_DEF  = 0.;
const ChemDriver::TRANSPORT transport_DEF = ChemDriver::CD_EG;
const std::string tranlib_linking_file="tran.asc";

namespace
{
    bool initialized = false;
    ChemDriver::TRANSPORT transport = transport_DEF;

    void ChemDriver_Finalize() {
      initialized = false;
      transport = transport_DEF;
    }
}

extern "C" {
  int FORT_USINGEG() {
    return (int)(transport == ChemDriver::CD_EG);
  }

  int FORT_USINGMC() {
    return (int)(transport == ChemDriver::CD_TRANLIB);
  }

  void FORT_MCLINKNAME(int* codedName, int* len, int* maxlen) {
    Array<int> file = ChemDriver::encodeStringForFortran(tranlib_linking_file);
    *len = file.size();
    if (*len > *maxlen) {
      BoxLib::Abort("Fort memory bust...need more space for tranlib linking file name");
    }

    for (int i=0; i<file.size(); ++i) {
      codedName[i] = file[i];
    }
  }

  struct ReactionData {
    double activation_units, prefactor_units, phase_units;
    double fwd_A,fwd_beta,fwd_Ea;
    double low_A,low_beta,low_Ea;
    double rev_A,rev_beta,rev_Ea;
    int is_PD, troe_len, sri_len;
    double troe_a,troe_Ts, troe_Tss, troe_Tsss;
    double sri_a, sri_b, sri_c, sri_d, sri_e;
  };

  struct ReactionData* GetReactionData(int id);
  struct ReactionData* GetDefaultReactionData(int id);
  void SetReactionData(int id, const struct ReactionData * rhs);
}

ChemDriver::ChemDriver ()
    :
    mHtoTerrMAX(HtoTerrMAX_DEF),
    mHtoTiterMAX(HtoTiterMAX_DEF)
{
    FORT_SETTMINTRANS(&Tmin_trans_DEF);

    //FORT_SETVERBOSEVODE();

    if (!initialized)
    {
        initOnce();
        BoxLib::ExecOnFinalize(ChemDriver_Finalize);
        initialized = true;
    }
    mTmpData.resize(mHtoTiterMAX);
}

ChemDriver::~ChemDriver ()
{
    FORT_FINALIZECHEM();
}

static 
Real* get_parameter_value(ChemDriver::REACTION_PARAMETER param_id,
                          ReactionData*                  r)
{
  // FIXME: Hard-coded list of parameters invites errors later...
  Real* ret = 0;
  if (param_id == ChemDriver::FWD_A)          {ret = &(r->fwd_A);}
  else if (param_id == ChemDriver::FWD_BETA)  {ret = &(r->fwd_beta);}
  else if (param_id == ChemDriver::FWD_EA)    {ret = &(r->fwd_Ea);}
  else if (param_id == ChemDriver::LOW_A)     {ret = &(r->low_A);}
  else if (param_id == ChemDriver::LOW_BETA)  {ret = &(r->low_beta);}
  else if (param_id == ChemDriver::LOW_EA)    {ret = &(r->low_Ea);}
  else if (param_id == ChemDriver::REV_A)     {ret = &(r->rev_A);}
  else if (param_id == ChemDriver::REV_BETA)  {ret = &(r->rev_beta);}
  else if (param_id == ChemDriver::REV_EA)    {ret = &(r->rev_Ea);}
  else if (param_id == ChemDriver::TROE_A)    {ret = &(r->troe_a);}
  else if (param_id == ChemDriver::TROE_TS)   {ret = &(r->troe_Ts);}
  else if (param_id == ChemDriver::TROE_TSS)  {ret = &(r->troe_Tss);}
  else if (param_id == ChemDriver::TROE_TSSS) {ret = &(r->troe_Tsss);}
  else if (param_id == ChemDriver::SRI_A)     {ret = &(r->sri_a);}
  else if (param_id == ChemDriver::SRI_B)     {ret = &(r->sri_b);}
  else if (param_id == ChemDriver::SRI_C)     {ret = &(r->sri_c);}
  else if (param_id == ChemDriver::SRI_D)     {ret = &(r->sri_d);}
  else if (param_id == ChemDriver::SRI_E)     {ret = &(r->sri_e);}
  else {
    BoxLib::Abort("Unknown reaction parameter");
  }
  return ret;
}

ChemDriver::Parameter::Parameter(int                _reaction_id,
                                 REACTION_PARAMETER _param_id)
{
  reaction_id = _reaction_id;
  param_id = _param_id;
}

Real
ChemDriver::Parameter::Parameter::Value() const
{
  return *get_parameter_value(param_id, GetReactionData(reaction_id));
}

Real
ChemDriver::Parameter::Parameter::DefaultValue() const
{
  return *get_parameter_value(param_id, GetDefaultReactionData(reaction_id));
}

void
ChemDriver::Parameter::Parameter::ResetToDefault()
{
  (*this) = *get_parameter_value(param_id, GetDefaultReactionData(reaction_id));
}

std::ostream&
ChemDriver::Parameter::operator<<(std::ostream& os) const
{
  // FIXME: Hard-coded list of parameters invites errors later...
  os << "(r: " << reaction_id << " p: ";
  if (param_id == FWD_A)          {os << "FWD_A";}
  else if (param_id == FWD_BETA)  {os << "FWD_BETA";}
  else if (param_id == FWD_EA)    {os << "FWD_EA";}
  else if (param_id == LOW_A)     {os << "LOW_A";}
  else if (param_id == LOW_BETA)  {os << "LOW_BETA";}
  else if (param_id == LOW_EA)    {os << "LOW_EA";}
  else if (param_id == REV_A)     {os << "REV_A";}
  else if (param_id == REV_BETA)  {os << "REV_BETA";}
  else if (param_id == REV_EA)    {os << "REV_EA";}
  else if (param_id == TROE_A)    {os << "TROE_A";}
  else if (param_id == TROE_TS)   {os << "TROE_TS";}
  else if (param_id == TROE_TSS)  {os << "TROE_TSS";}
  else if (param_id == TROE_TSSS) {os << "TROE_TSSS";}
  else if (param_id == SRI_A)     {os << "SRI_A";}
  else if (param_id == SRI_B)     {os << "SRI_B";}
  else if (param_id == SRI_C)     {os << "SRI_C";}
  else if (param_id == SRI_D)     {os << "SRI_D";}
  else if (param_id == SRI_E)     {os << "SRI_E";}
  else {
    BoxLib::Abort("Unknown reaction parameter");
  }
  os << " v: " << Value() << " DEF: " << DefaultValue() << ")";
  return os;
}

std::ostream& operator<< (std::ostream&  os, const ChemDriver::Parameter& param)
{
  return param.operator<<(os);
}

void
ChemDriver::Parameter::operator=(Real new_value)
{
  *get_parameter_value(param_id, GetReactionData(reaction_id)) = new_value;
}

void
ChemDriver::ResetAllParamsToDefault()
{
  // FIXME: Hard-coded list of parameters invites errors later...
  for (int i=0; i<numReactions(); ++i) {
    Parameter(i,FWD_A).ResetToDefault();
    Parameter(i,FWD_BETA).ResetToDefault();
    Parameter(i,FWD_EA).ResetToDefault();
    Parameter(i,LOW_A).ResetToDefault();
    Parameter(i,LOW_BETA).ResetToDefault();
    Parameter(i,LOW_EA).ResetToDefault();
    Parameter(i,REV_A).ResetToDefault();
    Parameter(i,REV_BETA).ResetToDefault();
    Parameter(i,REV_EA).ResetToDefault();
    Parameter(i,TROE_A).ResetToDefault();
    Parameter(i,TROE_TS).ResetToDefault();
    Parameter(i,TROE_TSS).ResetToDefault();
    Parameter(i,TROE_TSSS).ResetToDefault();
    Parameter(i,SRI_A).ResetToDefault();
    Parameter(i,SRI_B).ResetToDefault();
    Parameter(i,SRI_C).ResetToDefault();
    Parameter(i,SRI_D).ResetToDefault();
    Parameter(i,SRI_E).ResetToDefault();
  }
}


void
ChemDriver::SetTransport (const ChemDriver::TRANSPORT& tran_in)
{
  if (initialized) {
    BoxLib::Abort("Must set transport model before constructing the ChemDriver");
  }
  transport = tran_in;
}

ChemDriver::TRANSPORT
ChemDriver::Transport ()
{
  return transport;
}

void
ChemDriver::initOnce ()
{
    FORT_INITCHEM();
    getSpeciesNames();
    getElementNames();
    FORT_GETCKDIMPARAMS(&mMaxreac, &mMaxspec, &mMaxelts,  &mMaxord,
                        &mMaxthrdb, &mMaxtp,  &mMaxsp,    &mMaxspnml);
    getStoichCoeffs();

    ParmParse pp("ht");

    int  v_itol = 1;
    Real v_rtol = 1.e-10;
    Real v_atol = 1.e-10;

    pp.query("vode_itol",v_itol);
    pp.query("vode_rtol",v_rtol);
    pp.query("vode_atol",v_atol);

    BL_ASSERT(v_rtol > 0);
    BL_ASSERT(v_atol > 0);
    BL_ASSERT(v_itol == 1 || v_itol == 2);

    FORT_SETVODETOLS(&v_rtol,&v_atol,&v_itol);

    int  v_maxcyc = -1;

    pp.query("vode_max_subcycles",v_maxcyc);
    if (v_maxcyc > 0) {
      set_max_vode_subcycles(v_maxcyc);
    }
}

void
ChemDriver::set_verbose_vode()
{
    FORT_SETVERBOSEVODE();
}

void
ChemDriver::set_max_vode_subcycles(int maxcyc)
{
    FORT_SETVODESUBCYC(&maxcyc);
}

void
ChemDriver::set_species_Yscales(const std::string& scalesFile)
{
    Array<int> file = encodeStringForFortran(scalesFile);
    int len         = file.size();
    FORT_SETSPECSCALY(file.dataPtr(),&len);
}

void
ChemDriver::getSpeciesNames()
{
    int max_len = FORT_GETCKMAXNAMELEN();
    int* coded = new int[max_len];
    int Nspec = FORT_GETCKNUMSPEC();
    mSpeciesNames.clear();
    mSpeciesNames.resize(Nspec);
    for (int i=0; i<Nspec; ++i)
    {
        int ifort = i+1;
	int len = FORT_GETCKSPECNAME(&ifort, coded);
	mSpeciesNames[i] = decodeStringFromFortran(coded,len);
    }
    delete [] coded;
}

void
ChemDriver::getElementNames()
{
    const int max_len = 3;
    int* coded = new int[max_len];
    int Nelt = FORT_GETCKNUMELT();
    
    mElementNames.clear();
    mElementNames.resize(Nelt);
    for (int i=0; i<Nelt; ++i)
    {
        int ifort = i+1;
	int len = FORT_GETCKELTNAME(&ifort, coded);
	mElementNames[i] = decodeStringFromFortran(coded,len);
    }
    delete [] coded;
}

void
ChemDriver::getStoichCoeffs()
{
    mNu.resize(mMaxspec * mMaxreac);
    FORT_SETNU(mNu.dataPtr(),mNu.size());
}

Array<int>
ChemDriver::reactionsWithXonL(const std::string& specName) const
{
    const int idx = index(specName) + 1;
    int Nreacs = -1;
    Array<int> reactions(numReactions());
    FORT_FINDLHS(reactions.dataPtr(),&Nreacs,&idx);
    reactions.resize(Nreacs);
    for (int i=0; i<reactions.size(); ++i)
        reactions[i]--;
    return reactions;
}

Array<int>
ChemDriver::reactionsWithXonR(const std::string& specName) const
{
    const int idx = index(specName) + 1;
    int Nreacs = -1;
    Array<int> reactions(numReactions());
    FORT_FINDRHS(reactions.dataPtr(),&Nreacs,&idx);
    reactions.resize(Nreacs);
    for (int i=0; i<reactions.size(); ++i)
        reactions[i]--;
    return reactions;
}

int
ChemDriver::numberOfElementXinSpeciesY(const std::string& eltX,
                                          const std::string& spcY) const
{
    const int eltID = indexElt(eltX);
    const int spID  = index(spcY);
    return FORT_CKELTXINSPY(&eltID,&spID);
}

Array<std::pair<std::string,int> >
ChemDriver::specCoeffsInReactions(int reacIdx) const
{
    Array<int> KI(mMaxsp);
    Array<int> NU(mMaxsp);
    int Nids = 0;
    const int fortReacIdx = reacIdx + 1;
    FORT_CKINU(&Nids,KI.dataPtr(),&mMaxsp,NU.dataPtr(),&mMaxsp,&fortReacIdx,mNu.dataPtr());
    Array<std::pair<std::string,int> > result(Nids);
    for (int i=0; i<Nids; ++i)
        result[i] = std::pair<std::string,int>(mSpeciesNames[KI[i]-1],NU[i]);
    return result;
}

std::string
ChemDriver::reactionString(int reacIdx) const
{
    const int max_len = 72;
    const int fortReacIdx = reacIdx + 1;
    int* coded = new int[max_len];
    int len = FORT_CKSYMR(&fortReacIdx, coded);
    std::string line = decodeStringFromFortran(coded,len);
    delete [] coded;
    return line;
}

Array<Real>
ChemDriver::speciesMolecWt() const
{
    Array<Real> mwt(numSpecies());
    FORT_GETCKMWT(mwt.dataPtr());
    return mwt;
}

Array<Real>
ChemDriver::elementAtomicWt() const
{
    Array<Real> awt(numElements());
    FORT_GETCKAWT(awt.dataPtr());
    return awt;
}

Array<Real>
ChemDriver::massFracToMoleFrac(const Array<Real>& Y) const
{
    int nc = Y.size();
    BL_ASSERT(nc = numSpecies());
    Array<Real> X(nc);
    Box box(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(0,0,0)));
    FORT_MASSTOMOLE(box.loVect(), box.hiVect(),
		    Y.dataPtr(), ARLIM(box.loVect()), ARLIM(box.hiVect()),
		    X.dataPtr(), ARLIM(box.loVect()), ARLIM(box.hiVect()));
    return X;
}

Array<Real>
ChemDriver::moleFracToMassFrac(const Array<Real>& X) const
{
    int nc = X.size();
    BL_ASSERT(nc==numSpecies());
    Array<Real> Y(nc);
    Box box(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(0,0,0)));
    FORT_MOLETOMASS(box.loVect(), box.hiVect(),
		    X.dataPtr(), ARLIM(box.loVect()), ARLIM(box.hiVect()),
		    Y.dataPtr(), ARLIM(box.loVect()), ARLIM(box.hiVect()));
    return Y;
}

void
ChemDriver::normalizeMassFrac(FArrayBox&       Ynorm,
                                 const FArrayBox& Y,
                                 const std::string&   excessSpecies,
                                 const Box&       box,
                                 int              sCompY,
                                 int              sCompYn) const
{
    BL_ASSERT(sCompY+numSpecies() <= Y.nComp());
    BL_ASSERT(sCompYn+numSpecies() <= Ynorm.nComp());
    
    const Box& mabx = Y.box();
    const Box& mobx = Ynorm.box();
    
    Box ovlp = box & mabx & mobx;
    if (!ovlp.ok())
	return;

    int xsID = index(excessSpecies) + 1;
    BL_ASSERT(xsID > 0);
    FORT_NORMMASS(ovlp.loVect(), ovlp.hiVect(), &xsID,
                  Y.dataPtr(sCompY),      ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
                  Ynorm.dataPtr(sCompYn), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()));
}

int
ChemDriver::numReactions() const
{
    return FORT_GETCKNUMREAC();
}

void
ChemDriver::molarProduction(FArrayBox&       Q,
			       const std::string&   specName,
			       const FArrayBox& C,
			       const FArrayBox& T,
			       const Box&       box,
			       int              sCompC,
			       int              sCompT,
			       int              sCompQ) const
{
    const int Nreacs = numReactions();
    BL_ASSERT(Q.nComp() >= sCompQ + Nreacs);
    BL_ASSERT(C.nComp() >= sCompC + numSpecies());
    BL_ASSERT(T.nComp() > sCompT);

    const Box& mabx = C.box();
    const Box& mbbx = T.box();
    const Box& mobx = Q.box();
    const Box ovlp = box & mabx & mbbx & mobx;

    const int idx = index(specName) + 1; // to fortran indexing

    FORT_MOLPROD(ovlp.loVect(), ovlp.hiVect(), &idx,
		 Q.dataPtr(sCompQ), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()),
		 C.dataPtr(sCompC), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
		 T.dataPtr(sCompT), ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()));

}

void
ChemDriver::fwdRevReacRatesGivenXTP(FArrayBox&        FwdK,
                                       FArrayBox&        RevK,
                                       const Array<int>& rxnIDs,
                                       const FArrayBox&  X,
                                       const FArrayBox&  T,
                                       Real              Patm,
                                       const Box&        box,
                                       int               sCompX,
                                       int               sCompT,
                                       int               sCompFwdK,
                                       int               sCompRevK) const
{
    const int Nspec = numSpecies();
    const int Nreac = rxnIDs.size();
    BL_ASSERT(FwdK.nComp() >= sCompFwdK+Nreac);
    BL_ASSERT(RevK.nComp() >= sCompRevK+Nreac);
    BL_ASSERT(X.nComp() >= sCompX + Nspec);
    BL_ASSERT(T.nComp() > sCompT);
    
    const Box& mabx = X.box();
    const Box& mbbx = T.box();
    const Box& moabx = FwdK.box();
    const Box& mobbx = RevK.box();
    
    Box ovlp = box & mabx & mbbx & moabx & mobbx;
    if( ! ovlp.ok() ) return;
    
    FORT_FRrateXTP(ovlp.loVect(), ovlp.hiVect(),
                   X.dataPtr(sCompX),        ARLIM(mabx.loVect()),  ARLIM(mabx.hiVect()),
                   T.dataPtr(sCompT),        ARLIM(mbbx.loVect()),  ARLIM(mbbx.hiVect()),
                   FwdK.dataPtr(sCompFwdK), ARLIM(moabx.loVect()), ARLIM(moabx.hiVect()),
                   RevK.dataPtr(sCompRevK), ARLIM(mobbx.loVect()), ARLIM(mobbx.hiVect()),
                   &Patm,rxnIDs.dataPtr(),&Nreac);
}
    
void
ChemDriver::heatRelease(FArrayBox&       Q,
                           const FArrayBox& Y,
                           const FArrayBox& T,
                           Real             Patm,
                           const Box&       box,
                           int              sCompY,
                           int              sCompT,
                           int              sCompQ) const
{
    const int Nspec = numSpecies();
    BL_ASSERT(Q.nComp() > sCompQ);
    BL_ASSERT(Y.nComp() >= sCompY + Nspec);
    BL_ASSERT(T.nComp() > sCompT);

    const Box& mabx = Y.box();
    const Box& mbbx = T.box();
    const Box& mobx = Q.box();
    
    Box ovlp = box & mabx & mbbx & mobx;
    if( ! ovlp.ok() ) return;
    
    FORT_HTRLS(ovlp.loVect(), ovlp.hiVect(),
               Y.dataPtr(sCompY), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
               T.dataPtr(sCompT), ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
               Q.dataPtr(sCompQ), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()),
               &Patm);
}

void
ChemDriver::reactionRateY(FArrayBox&       Ydot,
                             const FArrayBox& Y,
                             const FArrayBox& T,
                             Real             Patm,
                             const Box&       box,
                             int              sCompY,
                             int              sCompT,
                             int              sCompYdot) const
{
    const int Nspec = numSpecies();
    BL_ASSERT(Ydot.nComp() >= sCompYdot + Nspec);
    BL_ASSERT(Y.nComp() >= sCompY + Nspec);
    BL_ASSERT(T.nComp() > sCompT);

    const Box& mabx = Y.box();
    const Box& mbbx = T.box();
    const Box& mobx = Ydot.box();
    
    Box ovlp = box & mabx & mbbx & mobx;
    if( ! ovlp.ok() ) return;
    
    FORT_RRATEY(ovlp.loVect(), ovlp.hiVect(),
                Y.dataPtr(sCompY),       ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
                T.dataPtr(sCompT),       ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
                Ydot.dataPtr(sCompYdot), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()),
                &Patm);
}

#ifdef LMC_SDC
void
ChemDriver::reactionRateRhoY(FArrayBox&       RhoYdot,
                             const FArrayBox& RhoY,
                             const FArrayBox& RhoH,
                             const FArrayBox& T,
                             Real             Patm,
                             const Box&       box,
                             int              sCompRhoY,
                             int              sCompRhoH,
                             int              sCompT,
                             int              sCompRhoYdot) const
{
    const int Nspec = numSpecies();
    BL_ASSERT(RhoYdot.nComp() >= sCompRhoYdot + Nspec);
    BL_ASSERT(RhoY.nComp() >= sCompRhoY + Nspec);
    BL_ASSERT(RhoH.nComp() > sCompRhoH);
    BL_ASSERT(T.nComp() > sCompT);

    const Box& mabx = RhoY.box();
    const Box& mbbx = RhoH.box();
    const Box& mcbx = T.box();
    const Box& mobx = RhoYdot.box();
    
    Box ovlp = box & mabx & mbbx & mcbx & mobx;
    if( ! ovlp.ok() ) return;
    
    FORT_RRATERHOY(ovlp.loVect(), ovlp.hiVect(),
                   RhoY.dataPtr(sCompRhoY),       ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
                   RhoH.dataPtr(sCompRhoH),       ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
                   T.dataPtr(sCompT),             ARLIM(mcbx.loVect()), ARLIM(mcbx.hiVect()),
                   RhoYdot.dataPtr(sCompRhoYdot), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()),
                   &Patm);
}
#endif

void
ChemDriver::massFracToMoleFrac(FArrayBox&       X,
				  const FArrayBox& Y,
				  const Box&       box,
				  int              sCompY,
				  int              sCompX) const
{
    BL_ASSERT(sCompY+numSpecies() <= Y.nComp());
    BL_ASSERT(sCompX+numSpecies() <= X.nComp());
    
    const Box& mabx = Y.box();
    const Box& mobx = X.box();
    
    Box ovlp = box & mabx & mobx;
    if (!ovlp.ok())
	return;
    
    FORT_MASSTOMOLE(ovlp.loVect(), ovlp.hiVect(),
		    Y.dataPtr(sCompY), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
		    X.dataPtr(sCompX), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()));
}

void
ChemDriver::moleFracToMassFrac(FArrayBox&       Y,
				  const FArrayBox& X,
				  const Box&       box,
				  int              sCompX,
				  int              sCompY) const
{
    BL_ASSERT(sCompX+numSpecies() <= X.nComp());
    BL_ASSERT(sCompY+numSpecies() <= Y.nComp());
    
    const Box& mobx = X.box();
    const Box& mabx = Y.box();
    
    Box ovlp = box & mobx & mabx;
    if (!ovlp.ok())
	return;
    
    FORT_MOLETOMASS(ovlp.loVect(), ovlp.hiVect(),
		    X.dataPtr(sCompX), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()),
		    Y.dataPtr(sCompY), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()));
}

void
ChemDriver::massFracToMolarConc(FArrayBox&       C,
				   const FArrayBox& Y,
				   const FArrayBox& T,
				   Real             Patm,
				   const Box&       box,
				   int              sCompY,
				   int              sCompT,
				   int              sCompC) const
{
    BL_ASSERT(sCompC+numSpecies() <= C.nComp());
    BL_ASSERT(sCompY+numSpecies() <= Y.nComp());
    BL_ASSERT(sCompT+1 <= T.nComp());
    
    const Box& mobx = C.box();
    const Box& mabx = Y.box();
    const Box& mbbx = T.box();
    
    Box ovlp = box & mabx & mbbx & mobx;
    if (!ovlp.ok())
	return;
    
    FORT_MASSTP_TO_CONC(ovlp.loVect(), ovlp.hiVect(), &Patm,
		    Y.dataPtr(sCompY), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
		    T.dataPtr(sCompT), ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
		    C.dataPtr(sCompC), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()));
}

void
ChemDriver::massFracToMolarConc(FArrayBox&       C,
				   const FArrayBox& Y,
				   const FArrayBox& T,
				   const FArrayBox& Rho,
				   const Box&       box,
				   int              sCompR,
				   int              sCompY,
				   int              sCompT,
				   int              sCompC) const
{
    BL_ASSERT(sCompC+numSpecies() <= C.nComp());
    BL_ASSERT(sCompY+numSpecies() <= Y.nComp());
    BL_ASSERT(sCompT+1 <= T.nComp());
    BL_ASSERT(sCompR <= Rho.nComp());

    const Box& mobx = C.box();
    const Box& mabx = Y.box();
    const Box& mbbx = T.box();
    const Box& mrbx = Rho.box();

    Box ovlp = box & mabx & mbbx & mobx & mrbx;
    if (!ovlp.ok())
	return;
    
    FORT_MASSR_TO_CONC(ovlp.loVect(), ovlp.hiVect(), 
		    Y.dataPtr(sCompY), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
		    T.dataPtr(sCompT), ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
		    Rho.dataPtr(sCompR),ARLIM(mrbx.loVect()),ARLIM(mrbx.hiVect()),
		    C.dataPtr(sCompC), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()));
}

void
ChemDriver::molarConcToMoleFrac(FArrayBox&       X,
                                   const FArrayBox& C,
                                   const Box&       box,
                                   int              sCompC,
                                   int              sCompX) const
{
    BL_ASSERT(sCompX+numSpecies() <= X.nComp());
    BL_ASSERT(sCompC+numSpecies() <= C.nComp());
    
    const Box& mobx = X.box();
    const Box& mabx = C.box();
    
    Box ovlp = box & mabx & mobx;
    if (!ovlp.ok())
	return;
    
    FORT_CONC_TO_MOLE(ovlp.loVect(), ovlp.hiVect(),
		      C.dataPtr(sCompC), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
		      X.dataPtr(sCompX), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()));
}

Array<int>
ChemDriver::encodeStringForFortran(const std::string& astr)
{
    long length = astr.size();
    Array<int> result(length);
    for (int i = 0; i < length; ++i)
        result[i] = astr[i];
    return result;
}

std::string
ChemDriver::decodeStringFromFortran(const int* coded,
				       int        length)
{
    std::string result;
    for (int i = 0; i < length; ++i)
        result += coded[i];
    return result;
}

#include "iostream"
using std::cout;
using std::endl;
void
ChemDriver::solveTransient(FArrayBox&        Ynew,
                           FArrayBox&        Tnew,
                           const FArrayBox&  Yold,
                           const FArrayBox&  Told,
                           FArrayBox&        FuncCount,
                           const Box&        box,
                           int               sCompY,
                           int               sCompT,
                           Real              dt,
                           Real              Patm,
                           FArrayBox*        chemDiag,
                           bool              use_stiff_solver) const
{
    BL_ASSERT(sCompY+numSpecies() <= Ynew.nComp());
    BL_ASSERT(sCompY+numSpecies() <= Yold.nComp());
    BL_ASSERT(sCompT < Tnew.nComp());
    BL_ASSERT(sCompT < Told.nComp());
    
    BL_ASSERT(Ynew.box().contains(box) && Yold.box().contains(box));
    BL_ASSERT(Tnew.box().contains(box) && Told.box().contains(box));

    const int do_diag  = (chemDiag!=0);
    Real*     diagData = do_diag ? chemDiag->dataPtr() : 0;
    const int do_stiff = (use_stiff_solver);
    FORT_CONPSOLV(box.loVect(), box.hiVect(),
                  Ynew.dataPtr(sCompY), ARLIM(Ynew.loVect()), ARLIM(Ynew.hiVect()),
                  Tnew.dataPtr(sCompT), ARLIM(Tnew.loVect()), ARLIM(Tnew.hiVect()),
                  Yold.dataPtr(sCompY), ARLIM(Yold.loVect()), ARLIM(Yold.hiVect()),
                  Told.dataPtr(sCompT), ARLIM(Told.loVect()), ARLIM(Told.hiVect()),
                  FuncCount.dataPtr(),
                  ARLIM(FuncCount.loVect()), ARLIM(FuncCount.hiVect()),
                  &Patm, &dt, diagData, &do_diag, &do_stiff);
}

#ifdef LMC_SDC
void
ChemDriver::solveTransient_sdc(FArrayBox&        rhoYnew,
			       FArrayBox&        rhoHnew,
			       FArrayBox&        Tnew,
			       const FArrayBox&  rhoYold,
			       const FArrayBox&  rhoHold,
			       const FArrayBox&  Told,
			       const FArrayBox&  const_src,
			       FArrayBox&        I_R,
			       FArrayBox&        FuncCount,
			       const Box&        box,
			       int               sComprhoY,
			       int               sComprhoH,
			       int               sCompT,
			       Real              dt,
			       Real              Patm,
			       FArrayBox*        chemDiag,
                               bool              use_stiff_solver) const
{
    BL_ASSERT(sComprhoY+numSpecies() <= rhoYnew.nComp());
    BL_ASSERT(sComprhoY+numSpecies() <= rhoYold.nComp());
    BL_ASSERT(sComprhoH < rhoHnew.nComp());
    BL_ASSERT(sComprhoH < rhoHold.nComp());
    BL_ASSERT(sCompT    < Tnew.nComp());
    BL_ASSERT(sCompT    < Told.nComp());
    
    BL_ASSERT(rhoYnew.box().contains(box) && rhoYold.box().contains(box));
    BL_ASSERT(rhoHnew.box().contains(box) && rhoHold.box().contains(box));
    BL_ASSERT(Tnew.box().contains(box)    && Told.box().contains(box));

    const int do_diag  = (chemDiag!=0);
    Real*     diagData = do_diag ? chemDiag->dataPtr() : 0;
    const int do_stiff = (use_stiff_solver);

    FORT_CONPSOLV_SDC(box.loVect(), box.hiVect(),
                      rhoYnew.dataPtr(sComprhoY), ARLIM(rhoYnew.loVect()),   ARLIM(rhoYnew.hiVect()),
                      rhoHnew.dataPtr(sComprhoH), ARLIM(rhoHnew.loVect()),   ARLIM(rhoHnew.hiVect()),
                      Tnew.dataPtr(sCompT),       ARLIM(Tnew.loVect()),      ARLIM(Tnew.hiVect()),
                      rhoYold.dataPtr(sComprhoY), ARLIM(rhoYold.loVect()),   ARLIM(rhoYold.hiVect()),
                      rhoHold.dataPtr(sComprhoH), ARLIM(rhoHold.loVect()),   ARLIM(rhoHold.hiVect()),
                      Told.dataPtr(sCompT),       ARLIM(Told.loVect()),      ARLIM(Told.hiVect()),
                      const_src.dataPtr(0),       ARLIM(const_src.loVect()), ARLIM(const_src.hiVect()),
		      I_R.dataPtr(0),             ARLIM(I_R.loVect()),       ARLIM(I_R.hiVect()),
                      FuncCount.dataPtr(),        ARLIM(FuncCount.loVect()), ARLIM(FuncCount.hiVect()),
		      &Patm, &dt, diagData, &do_diag, &do_stiff);
}
#endif

void
ChemDriver::getMixAveragedRhoDiff(FArrayBox&       rhoD,
                                     const FArrayBox& Y,
                                     const FArrayBox& T,
                                     Real             Patm,
                                     const Box&       box,
                                     int              sCompY,
                                     int              sCompT,
                                     int              sCompRD) const
{
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    BL_ASSERT(rhoD.nComp() >= sCompRD+numSpecies());
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(rhoD.box().contains(box));

    int do_temp    = 0;
    int do_VelVisc = 0;

    FORT_MIXAVG_RHODIFF_TEMP(box.loVect(), box.hiVect(),
                             rhoD.dataPtr(sCompRD),ARLIM(rhoD.loVect()),ARLIM(rhoD.hiVect()),
                             T.dataPtr(sCompT),    ARLIM(T.loVect()),   ARLIM(T.hiVect()),
                             Y.dataPtr(sCompY),    ARLIM(Y.loVect()),   ARLIM(Y.hiVect()),
                             &Patm,&do_temp,&do_VelVisc);


}

void
ChemDriver::getMixShearVisc(FArrayBox&       eta,
                               const FArrayBox& Y,
                               const FArrayBox& T,
                               const Box&       box,
                               int              sCompY,
                               int              sCompT,
                               int              sCompE) const
{
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    BL_ASSERT(eta.nComp() >= sCompE);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(eta.box().contains(box));
    
    FORT_MIX_SHEAR_VISC(box.loVect(), box.hiVect(),
                        eta.dataPtr(sCompE),ARLIM(eta.loVect()),ARLIM(eta.hiVect()),
                        T.dataPtr(sCompT),  ARLIM(T.loVect()),  ARLIM(T.hiVect()),
                        Y.dataPtr(sCompY),  ARLIM(Y.loVect()),  ARLIM(Y.hiVect()));
}

void
ChemDriver::getRhoGivenPTY(FArrayBox&       Rho,
			      Real             Patm,
			      const FArrayBox& T,
			      const FArrayBox& Y,
			      const Box&       box,
			      int              sCompT,
			      int              sCompY,
			      int              sCompR) const
{
    BL_ASSERT(Rho.nComp() > sCompR);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());

    BL_ASSERT(Rho.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_RHOfromPTY(box.loVect(), box.hiVect(),
		    Rho.dataPtr(sCompR), ARLIM(Rho.loVect()), ARLIM(Rho.hiVect()),
		    T.dataPtr(sCompT),   ARLIM(T.loVect()),   ARLIM(T.hiVect()),
		    Y.dataPtr(sCompY),   ARLIM(Y.loVect()),   ARLIM(Y.hiVect()),
		    &Patm);
}

void
ChemDriver::getRhoGivenPvTY(FArrayBox&       Rho,
			      const FArrayBox& P,
			      const FArrayBox& T,
			      const FArrayBox& Y,
			      const Box&       box,
			      int              sCompP,
			      int              sCompT,
			      int              sCompY,
			      int              sCompR) const
{
    BL_ASSERT(Rho.nComp() > sCompR);
    BL_ASSERT(P.nComp() > sCompP);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());

    BL_ASSERT(Rho.box().contains(box));
    BL_ASSERT(P.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_RHOfromPvTY(box.loVect(), box.hiVect(),
		    Rho.dataPtr(sCompR), ARLIM(Rho.loVect()), ARLIM(Rho.hiVect()),
		    T.dataPtr(sCompT),   ARLIM(T.loVect()),   ARLIM(T.hiVect()),
		    Y.dataPtr(sCompY),   ARLIM(Y.loVect()),   ARLIM(Y.hiVect()),
		    P.dataPtr(sCompP),   ARLIM(P.loVect()),   ARLIM(P.hiVect()));
}

void
ChemDriver::getPGivenRTY(FArrayBox&       p,
			    const FArrayBox& Rho,
			    const FArrayBox& T,
			    const FArrayBox& Y,
			    const Box&       box,
			    int              sCompR,
			    int              sCompT,
			    int              sCompY,
			    int              sCompP) const
{
    BL_ASSERT(p.nComp() > sCompP);
    BL_ASSERT(Rho.nComp() > sCompR);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    
    BL_ASSERT(p.box().contains(box));
    BL_ASSERT(Rho.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_PfromRTY(box.loVect(), box.hiVect(),
		  p.dataPtr(sCompP),   ARLIM(p.loVect()),   ARLIM(p.hiVect()),
		  Rho.dataPtr(sCompR), ARLIM(Rho.loVect()), ARLIM(Rho.hiVect()),
		  T.dataPtr(sCompT),   ARLIM(T.loVect()),   ARLIM(T.hiVect()),
		  Y.dataPtr(sCompY),   ARLIM(Y.loVect()),   ARLIM(Y.hiVect()));
}

void
ChemDriver::getTGivenPRY(FArrayBox&       T,
                         Real             Patm,
                         const FArrayBox& Rho,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompR,
                         int              sCompY,
                         int              sCompT) const
{
    BL_ASSERT(Rho.nComp() > sCompR);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    
    BL_ASSERT(Rho.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_TfromPRY(box.loVect(), box.hiVect(),
		  T.dataPtr(sCompT),   ARLIM(T.loVect()),   ARLIM(T.hiVect()),
		  Rho.dataPtr(sCompR), ARLIM(Rho.loVect()), ARLIM(Rho.hiVect()),
		  Y.dataPtr(sCompY),   ARLIM(Y.loVect()),   ARLIM(Y.hiVect()),
                  &Patm);
}

void
ChemDriver::getCpmixGivenTY(FArrayBox&       cpmix,
			       const FArrayBox& T,
			       const FArrayBox& Y,
			       const Box&       box,
			       int              sCompT,
			       int              sCompY,
			       int              sCompCp) const
{
    BL_ASSERT(cpmix.nComp() > sCompCp);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());

    BL_ASSERT(cpmix.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_CPMIXfromTY(box.loVect(), box.hiVect(),
		     cpmix.dataPtr(sCompCp),ARLIM(cpmix.loVect()), ARLIM(cpmix.hiVect()),
		     T.dataPtr(sCompT),     ARLIM(T.loVect()),     ARLIM(T.hiVect()),
		     Y.dataPtr(sCompY),     ARLIM(Y.loVect()),     ARLIM(Y.hiVect()));
}

void
ChemDriver::getCvmixGivenTY(FArrayBox&       cvmix,
			       const FArrayBox& T,
			       const FArrayBox& Y,
			       const Box&       box,
			       int              sCompT,
			       int              sCompY,
			       int              sCompCv) const
{
    BL_ASSERT(cvmix.nComp() > sCompCv);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());

    BL_ASSERT(cvmix.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_CVMIXfromTY(box.loVect(), box.hiVect(),
		     cvmix.dataPtr(sCompCv),ARLIM(cvmix.loVect()), ARLIM(cvmix.hiVect()),
		     T.dataPtr(sCompT),     ARLIM(T.loVect()),     ARLIM(T.hiVect()),
		     Y.dataPtr(sCompY),     ARLIM(Y.loVect()),     ARLIM(Y.hiVect()));
}

void
ChemDriver::getHmixGivenTY(FArrayBox&       hmix,
			      const FArrayBox& T,
			      const FArrayBox& Y,
			      const Box&       box,
			      int              sCompT,
			      int              sCompY,
			      int              sCompH) const
{
    BL_ASSERT(hmix.nComp() > sCompH);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    
    BL_ASSERT(hmix.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_HMIXfromTY(box.loVect(), box.hiVect(),
		    hmix.dataPtr(sCompH),ARLIM(hmix.loVect()), ARLIM(hmix.hiVect()),
		    T.dataPtr(sCompT),   ARLIM(T.loVect()),    ARLIM(T.hiVect()),
		    Y.dataPtr(sCompY),   ARLIM(Y.loVect()),    ARLIM(Y.hiVect()));
}

void
ChemDriver::getMwmixGivenY(FArrayBox&       mwmix,
			      const FArrayBox& Y,
			      const Box&       box,
			      int              sCompY,
			      int              sCompMw) const
{
    BL_ASSERT(mwmix.nComp() > sCompMw);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());

    BL_ASSERT(mwmix.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_MWMIXfromY(box.loVect(), box.hiVect(),
		    mwmix.dataPtr(sCompMw),ARLIM(mwmix.loVect()), ARLIM(mwmix.hiVect()),
		    Y.dataPtr(sCompY),     ARLIM(Y.loVect()),     ARLIM(Y.hiVect()));
}

void
ChemDriver::getCpGivenT(FArrayBox&       cp,
			   const FArrayBox& T,
			   const Box&       box,
			   int              sCompT,
			   int              sCompCp) const
{
    BL_ASSERT(cp.nComp() >= sCompCp + numSpecies());
    BL_ASSERT(T.nComp() > sCompT);

    BL_ASSERT(cp.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    
    FORT_CPfromT(box.loVect(), box.hiVect(),
		 cp.dataPtr(sCompCp), ARLIM(cp.loVect()), ARLIM(cp.hiVect()),
		 T.dataPtr(sCompT),   ARLIM(T.loVect()),  ARLIM(T.hiVect()));
}

void
ChemDriver::getHGivenT(FArrayBox&       h,
			  const FArrayBox& T,
			  const Box&       box,
			  int              sCompT,
			  int              sCompH) const
{
    BL_ASSERT(h.nComp() >= sCompH + numSpecies());
    BL_ASSERT(T.nComp() > sCompT);

    BL_ASSERT(h.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    
    FORT_HfromT(box.loVect(), box.hiVect(),
		h.dataPtr(sCompH), ARLIM(h.loVect()), ARLIM(h.hiVect()),
		T.dataPtr(sCompT), ARLIM(T.loVect()), ARLIM(T.hiVect()));
}

int
ChemDriver::getTGivenHY(FArrayBox&       T,
                        const FArrayBox& H,
                        const FArrayBox& Y,
                        const Box&       box,
                        int              sCompH,
                        int              sCompY,
                        int              sCompT,
                        const Real&      errMAX) const
{
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(H.nComp() > sCompH);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(H.box().contains(box));
    BL_ASSERT(Y.box().contains(box));

    Real solveTOL = (errMAX<0 ? mHtoTerrMAX : errMAX);
    
    return FORT_TfromHY(box.loVect(), box.hiVect(),
			T.dataPtr(sCompT), ARLIM(T.loVect()), ARLIM(T.hiVect()),
			H.dataPtr(sCompH), ARLIM(H.loVect()), ARLIM(H.hiVect()),
			Y.dataPtr(sCompY), ARLIM(Y.loVect()), ARLIM(Y.hiVect()),
			&solveTOL,&mHtoTiterMAX,mTmpData.dataPtr());
}

void
ChemDriver::getElementMoles(FArrayBox&       C_elt,
                               const std::string&   name,
                               const FArrayBox& C,
                               const Box&       box,
                               int              sCompC,
                               int              sCompC_elt) const
{
    BL_ASSERT(C.nComp() >= sCompC+numSpecies());
    BL_ASSERT(C_elt.nComp() > sCompC_elt);

    Array<int> name_enc = encodeStringForFortran(name);
    const int name_len = name_enc.size();

    FORT_GETELTMOLES(name_enc.dataPtr(), &name_len,
                     box.loVect(), box.hiVect(),
                     C_elt.dataPtr(sCompC_elt),
                     ARLIM(C_elt.loVect()),ARLIM(C_elt.hiVect()),
                     C.dataPtr(sCompC),ARLIM(C.loVect()),ARLIM(C.hiVect()));
}

Real
ChemDriver::getRuniversal() const
{
    return FORT_RUNIV();
}

Real
ChemDriver::getP1atm_MKS() const
{
  return FORT_P1ATMMKS();
}

void
ChemDriver::getOTradLoss_TDF(FArrayBox&       Qloss,
                                const FArrayBox& T,
                                const FArrayBox& X,
                                const Real       Patm,
                                const Real       T_bg,
                                const Box&       box,
                                int              sCompX,
                                int              sCompT,
                                int              sCompQ) const
{
    BL_ASSERT(Qloss.nComp() > sCompQ);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(X.nComp() >= sCompX + numSpecies());

    BL_ASSERT(Qloss.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(X.box().contains(box));

    FORT_OTrad_TDF(box.loVect(), box.hiVect(),
                   Qloss.dataPtr(sCompQ), ARLIM(Qloss.loVect()), ARLIM(Qloss.hiVect()),
                   T.dataPtr(sCompT),     ARLIM(T.loVect()),     ARLIM(T.hiVect()),
                   X.dataPtr(sCompX),     ARLIM(X.loVect()),     ARLIM(X.hiVect()),
                   &Patm, &T_bg);
}

ChemDriver::Edge::Edge (const std::string& n1, const std::string& n2, const Array<std::pair<int,Real> > rwl)
    : sp1(n1), sp2(n2), RWL(rwl) {}

ChemDriver::Edge::Edge (const std::string& n1, const std::string& n2, int reac, Real weight )
    : sp1(n1), sp2(n2) { RWL.push_back(std::pair<int,Real>(reac,weight)); }

int
ChemDriver::Edge::equivSign (const Edge& rhs) const
{   
    if ( (sp1 == rhs.sp1) && (sp2 == rhs.sp2) )
        return 1;
    
    else
        if ( (sp1 == rhs.sp2) && (sp2 == rhs.sp1) )
            return -1;
    
    return 0;
}

void
ChemDriver::Edge::combine (const Edge& rhs, const int   sgn)
{
    if (sgn!=0)
    {            
        int oldSize = RWL.size();
        int delSize = rhs.RWL.size();
        RWL.resize(oldSize + delSize);
        for (int i=0; i<delSize; ++i)
            RWL[oldSize+i] = std::pair<int,Real>(rhs.RWL[i].first,sgn*rhs.RWL[i].second);
    }
}

bool
ChemDriver::Edge::touchesSp(const std::string& rhs) const
{
    return ( (sp1 == rhs) || (sp2 == rhs) );
}

void
ChemDriver::Edge::reverse()
{
    for (int i=0; i<RWL.size(); ++i)
        RWL[i].second = - RWL[i].second;
    std::swap(sp1,sp2);
}

const Array<std::pair<int,Real> >&
ChemDriver::Edge::rwl () const
{
    return RWL;
}

const std::string
ChemDriver::Edge::left() const
{
    return sp1;
}

const std::string
ChemDriver::Edge::right() const
{
    return sp2;
}

bool
ChemDriver::Edge::operator<(const Edge& rhs) const
{
    return (this->left() == rhs.left() ?  this->right() < rhs.right() : this->left() < rhs.left());
}

std::ostream& operator<< (std::ostream& os, const ChemDriver::Edge& e)
{
    os << e.sp1 << " <=> " << e.sp2 << "  ( ";
    for (Array<std::pair<int,Real> >::const_iterator it=e.RWL.begin(); it!=e.RWL.end(); ++it)
        os << it->first << ":" << it->second << " ";
    os << ")";
    return os;
}

ChemDriver::Group::Group (const std::map<std::string,int>& eltCnts)
{
    mEltCnts = eltCnts;
}

ChemDriver::Group::Group (const ChemDriver::Group& rhs)
{
    mEltCnts = rhs.mEltCnts;
}

ChemDriver::Group
ChemDriver::Group::operator- (const ChemDriver::Group& rhs) const
{
    ChemDriver::Group val(*this);
    for (std::map<std::string,int>::const_iterator rit = rhs.mEltCnts.begin();rit!=rhs.mEltCnts.end();++rit)
    {
        std::map<std::string,int>::iterator fit = val.mEltCnts.find( rit->first );
        if ( fit != val.mEltCnts.end() )
        {
            int res = fit->second - rit->second;
            if (res != 0)
                fit->second = res;
            else
                val.mEltCnts.erase(fit);
        }
        else
            val.mEltCnts[rit->first] = -rit->second;
    }
    return val;
}

ChemDriver::Group
ChemDriver::Group::operator* (int rhs) const
{
    ChemDriver::Group val(*this);
    for (std::map<std::string,int>::iterator vit = val.mEltCnts.begin();vit!=val.mEltCnts.end();++vit)
        vit->second = rhs*vit->second;
    return val;
}

int
ChemDriver::Group::operator[](const std::string& id) const
{
    std::map<std::string,int>::const_iterator fit = mEltCnts.find( id );
    if (fit!=mEltCnts.end())
        return fit->second;
    return 0;
}

bool
ChemDriver::Group::sameSign () const
{
    if (mEltCnts.size()>0)
    {
        std::map<std::string,int>::const_iterator it = mEltCnts.begin();
        if (it->second < 0) {
            for (++it; it!=mEltCnts.end(); ++it)
                if (it->second > 0)
                    return false;
        }
        else
        {
            for (++it; it!=mEltCnts.end(); ++it)
                if (it->second < 0)
                    return false;
        }
    }
    return true;
}

bool
ChemDriver::Group::contains (const std::string& id) const
{
    return mEltCnts.find(id)!=mEltCnts.end();
}
    
Real
ChemDriver::Group::awt ()
{
    if (AtomicWeight.size() == 0)
        FillAtomicWeights(); // read up static atomic weight data
    Real myWt = 0;
    for (std::map<std::string,int>::const_iterator it = mEltCnts.begin(); it!=mEltCnts.end(); ++it)
    {
        Real thisAwt = AtomicWeight[it->first];
        myWt += std::abs(it->second) * thisAwt;
    }
    return myWt;
}

int
ChemDriver::Group::size() const
{
    int mySize = 0;
    for (std::map<std::string,int>::const_iterator it = mEltCnts.begin(); it!=mEltCnts.end(); ++it)
        mySize += std::abs(it->second);
    return mySize;
}

ChemDriver::Group operator* (int n, const ChemDriver::Group& g)
{
    return g.operator*(n);
}

std::map<std::string,Real> ChemDriver::Group::AtomicWeight; // Static initializer

void
ChemDriver::Group::FillAtomicWeights ()
{
    AtomicWeight["H"] = 1.00797;
    AtomicWeight["HE"] = 4.00260;
    AtomicWeight["LI"] = 6.93900;
    AtomicWeight["BE"] = 9.01220;
    AtomicWeight["B"] = 10.81100;
    AtomicWeight["C"] = 12.01115;
    AtomicWeight["N"] = 14.00670;
    AtomicWeight["O"] = 15.99940;
    AtomicWeight["F"] = 18.99840;
    AtomicWeight["NE"] = 20.18300;
    AtomicWeight["NA"] = 22.98980;
    AtomicWeight["MG"] = 24.31200;
    AtomicWeight["AL"] = 26.98150;
    AtomicWeight["SI"] = 28.08600;
    AtomicWeight["P"] = 30.97380;
    AtomicWeight["S"] = 32.06400;
    AtomicWeight["CL"] = 35.45300;
    AtomicWeight["AR"] = 39.94800;
    AtomicWeight["K"] = 39.10200;
    AtomicWeight["CA"] = 40.08000;
    AtomicWeight["SC"] = 44.95600;
    AtomicWeight["TI"] = 47.90000;
    AtomicWeight["V"] = 50.94200;
    AtomicWeight["CR"] = 51.99600;
    AtomicWeight["MN"] = 54.93800;
    AtomicWeight["FE"] = 55.84700;
    AtomicWeight["CO"] = 58.93320;
    AtomicWeight["NI"] = 58.71000;
    AtomicWeight["CU"] = 63.54000;
    AtomicWeight["ZN"] = 65.37000;
    AtomicWeight["GA"] = 69.72000;
    AtomicWeight["GE"] = 72.59000;
    AtomicWeight["AS"] = 74.92160;
    AtomicWeight["SE"] = 78.96000;
    AtomicWeight["BR"] = 79.90090;
    AtomicWeight["KR"] = 83.80000;
    AtomicWeight["RB"] = 85.47000;
    AtomicWeight["SR"] = 87.62000;
    AtomicWeight["Y"] = 88.90500;
    AtomicWeight["ZR"] = 91.22000;
    AtomicWeight["NB"] = 92.90600;
    AtomicWeight["MO"] = 95.94000;
    AtomicWeight["TC"] = 99.00000;
    AtomicWeight["RU"] = 101.07000;
    AtomicWeight["RH"] = 102.90500;
    AtomicWeight["PD"] = 106.40000;
    AtomicWeight["AG"] = 107.87000;
    AtomicWeight["CD"] = 112.40000;
    AtomicWeight["IN"] = 114.82000;
    AtomicWeight["SN"] = 118.69000;
    AtomicWeight["SB"] = 121.75000;
    AtomicWeight["TE"] = 127.60000;
    AtomicWeight["I"] = 126.90440;
    AtomicWeight["XE"] = 131.30000;
    AtomicWeight["CS"] = 132.90500;
    AtomicWeight["BA"] = 137.34000;
    AtomicWeight["LA"] = 138.91000;
    AtomicWeight["CE"] = 140.12000;
    AtomicWeight["PR"] = 140.90700;
    AtomicWeight["ND"] = 144.24000;
    AtomicWeight["PM"] = 145.00000;
    AtomicWeight["SM"] = 150.35000;
    AtomicWeight["EU"] = 151.96000;
    AtomicWeight["GD"] = 157.25000;
    AtomicWeight["TB"] = 158.92400;
    AtomicWeight["DY"] = 162.50000;
    AtomicWeight["HO"] = 164.93000;
    AtomicWeight["ER"] = 167.26000;
    AtomicWeight["TM"] = 168.93400;
    AtomicWeight["YB"] = 173.04000;
    AtomicWeight["LU"] = 174.99700;
    AtomicWeight["HF"] = 178.49000;
    AtomicWeight["TA"] = 180.94800;
    AtomicWeight["W"] = 183.85000;
    AtomicWeight["RE"] = 186.20000;
    AtomicWeight["OS"] = 190.20000;
    AtomicWeight["IR"] = 192.20000;
    AtomicWeight["PT"] = 195.09000;
    AtomicWeight["AU"] = 196.96700;
    AtomicWeight["HG"] = 200.59000;
    AtomicWeight["TL"] = 204.37000;
    AtomicWeight["PB"] = 207.19000;
    AtomicWeight["BI"] = 208.98000;
    AtomicWeight["PO"] = 210.00000;
    AtomicWeight["AT"] = 210.00000;
    AtomicWeight["RN"] = 222.00000;
    AtomicWeight["FR"] = 223.00000;
    AtomicWeight["RA"] = 226.00000;
    AtomicWeight["AC"] = 227.00000;
    AtomicWeight["TH"] = 232.03800;
    AtomicWeight["PA"] = 231.00000;
    AtomicWeight["U"] = 238.03000;
    AtomicWeight["NP"] = 237.00000;
    AtomicWeight["PU"] = 242.00000;
    AtomicWeight["AM"] = 243.00000;
    AtomicWeight["CM"] = 247.00000;
    AtomicWeight["BK"] = 249.00000;
    AtomicWeight["CF"] = 251.00000;
    AtomicWeight["ES"] = 254.00000;
    AtomicWeight["FM"] = 253.00000;
    AtomicWeight["D"] = 002.01410;
    AtomicWeight["E"] = 5.48578E-4;
}

std::ostream& operator<< (std::ostream& os, const ChemDriver::Group& g)
{
    os << "Group < ";
    for (std::map<std::string,int>::const_iterator it=g.mEltCnts.begin(); it!=g.mEltCnts.end(); ++it)
        os << it->first << ":" << it->second << " ";
    os << ">";
    return os;
}

std::list<ChemDriver::Edge>
ChemDriver::getEdges (const std::string& trElt, int PrintVerbose, int HackSplitting) const
{
    std::list<Edge> edges;
    std::map<std::string,Group> groups;
    const Array<std::string>& spNames = speciesNames();
    const Array<std::string>& eltNames = elementNames();
    for (int i=0; i<spNames.size(); ++i)
    {
        const std::string sp = spNames[i];
        std::map<std::string,int> eltCnts;
        for (int j=0; j<eltNames.size(); ++j)
        {
            const std::string elt = eltNames[j];
            int m = numberOfElementXinSpeciesY(elt,sp);
            if (m>0)
                eltCnts[elt] = m;

            groups[sp] = Group(eltCnts);
        }
    }
    for (int r=0; r<numReactions(); ++r)
    {
        const Array<std::pair<std::string,int> >& coeffs = specCoeffsInReactions(r);

        std::map<std::string,int> net;
        for (int i=0; i<coeffs.size(); ++i)
        {
            const std::string& sp=coeffs[i].first;
            const int co=coeffs[i].second;
            std::map<std::string,int>::iterator it = net.find(sp);
            if (it!=net.end())
            {
                net[sp] += co;
                if (net[sp]==0)
                    net.erase(it);
            }
            else
            {
                net[sp] = co;
            }
        }

        std::map<std::string,int> reac, prod;

        for (std::map<std::string,int>::const_iterator it = net.begin(); it!=net.end(); ++it)
        {
            const std::string& sp = it->first;
            const int n = it->second;
            if (n < 0)
            {
                if (numberOfElementXinSpeciesY(trElt,sp)>0)
                    reac[sp] = -n;
            }
            else
            {
                if (numberOfElementXinSpeciesY(trElt,sp)>0)
                    prod[sp] = n;
            }
        }

        int LR = reac.size();
        int LP = prod.size();
        
        if (LR==0 || LP==0)  // no tr-containing species in reac
            continue;
        
        if ((LR == 1) || (LP == 1)) // exactly one tr-containing species either side
        {
            for (std::map<std::string,int>::const_iterator rit=reac.begin(); rit!=reac.end(); ++rit)
            {
                const std::string& spcr = rit->first;
                const int cor = rit->second;

                for (std::map<std::string,int>::const_iterator pit=prod.begin(); pit!=prod.end(); ++pit)
                {
                    const std::string& spcp = pit->first;
                    const int cop = pit->second;
                    int w = std::min(cor*groups[spcr][trElt],cop*groups[spcp][trElt]);
                    edges.push_back(Edge(spcr,spcp,r,w));
                }
            }
            continue;
        }

        if ( (LR==2) && (LP==2) )  // two tr-containing species each side
        {
            std::map<std::string,int>::const_iterator r0 = reac.begin();
            std::map<std::string,int>::const_iterator r1 = r0; r1++;
            std::map<std::string,int>::const_iterator p0 = prod.begin();
            std::map<std::string,int>::const_iterator p1 = p0; p1++;
            
            const std::string& rs0 = r0->first;
            const std::string& rs1 = r1->first;
            const std::string& ps0 = p0->first;
            const std::string& ps1 = p1->first;
            
            int rc0 = r0->second;
            int rc1 = r1->second;
            int pc0 = p0->second;
            int pc1 = p1->second;

            Group b0(pc0 * groups[ps0] - rc0 * groups[rs0]);
            Group b1(pc1 * groups[ps1] - rc0 * groups[rs0]);
            int pick = 0;

            // HACK
            if ((HackSplitting==1) && (trElt=="H") && (r==61))
            {
                pick = 0;
                if (PrintVerbose>0)
                    cout << "...resorting to HACK!!" << endl;
            }
            else
            {
                if ((HackSplitting==1) && (trElt=="H") && (r==310))
                {
                    pick = 0;
                    if (PrintVerbose>0)
                        cout << "...resorting to HACK!!" << endl;
                }
                else
                {
                    if (b0.sameSign() && b1.sameSign())
                    {
                        if (b1.size() < b0.size())
                        {
                            pick = 1;
                            if (PrintVerbose>0)
                                cout << "...resorting to atom cnt" << endl;
                        }
                        if (b1.size() == b0.size())
                        {
                            if (b0.awt() > b1.awt())
                            {
                                pick = 1;
                                if (PrintVerbose>0)
                                    cout << "...resorting to tot awt" << endl;
                            }
                            else if ( (PrintVerbose>0) && (b0.awt() < b1.awt()) )
                            {
                                cout << "...resorting to tot awt" << endl;
                            }
                        }
                    }
                    else
                    {
                        if (b1.sameSign())
                            pick = 1;
                    }
                }
            }
            
            int nR0 = rc0*groups[rs0][trElt];
            int nR1 = rc1*groups[rs1][trElt];
            int nP0 = pc0*groups[ps0][trElt];
            int nP1 = pc1*groups[ps1][trElt];
            
            if (pick == 0)
            {
                edges.push_back(Edge(rs0,ps0,r,std::min(nR0,nP0)));
                if (nP0<nR0)
                    edges.push_back(Edge(rs0,ps1,r,nR0-nP0));
                
                edges.push_back(Edge(rs1,ps1,r,std::min(nR1,nP1)));
                if (nR0<nP0)
                    edges.push_back(Edge(rs1,ps0,r,nP0-nR0));
            }
            else
            {
                edges.push_back(Edge(rs0,ps1,r,std::min(nR0,nP1)));
                if (nP1<nR0)
                    edges.push_back(Edge(rs0,ps0,r,nR0-nP1));
                
                edges.push_back(Edge(rs1,ps0,r,std::min(nR1,nP0)));
                if (nR0<nP1)
                    edges.push_back(Edge(rs1,ps1,r,nP1-nR0));
            }
            continue;
        }

        cout << "Cannot decompose rxn: " << r << " " << LR << " " << LP << endl;
    }

    // Uniquify/combine edges
    std::list<Edge> uEdges;
    while (edges.size()!=0)
    {
        const Edge oe(edges.front()); edges.pop_front();
        bool foundIt = false;
        // There's probably a faster way to search this list...
        for (std::list<Edge>::iterator it=uEdges.begin(); !foundIt && it!=uEdges.end(); ++it)
        {
            int sgn = it->equivSign(oe);
            if (sgn!=0)
            {
                it->combine(oe,sgn);
                foundIt = true;
            }
        }
        if (!foundIt)
            uEdges.push_back(oe);
    }
    
    return uEdges;
}
