#include <winstd.H>

#include "ChemDriver.H"
#include "ChemDriver_F.H"
#include <ParallelDescriptor.H>

bool ChemDriver::mInitialized = false;

const Real HtoTerrMAX_DEF = 1.e-8;
const int HtoTiterMAX_DEF = 20;
const Real Tmin_trans_DEF = 0.;

ChemDriver::ChemDriver (const std::string& TransportFile)
    :
    mHtoTerrMAX(HtoTerrMAX_DEF),
    mHtoTiterMAX(HtoTiterMAX_DEF)
{
    FORT_SETTMINTRANS(&Tmin_trans_DEF);

    //FORT_SETVERBOSEVODE();

    if (!mInitialized)
    {
        initOnce(TransportFile);
        mInitialized = true;
    }
    mTmpData.resize(mHtoTiterMAX);
}

void
ChemDriver::initOnce (const std::string& TransportFile)
{
    const int NPROCS = ParallelDescriptor::NProcs();
    const int MYPROC = ParallelDescriptor::MyProc();
    //
    // Serialize access to files.
    // Some systems have severe limits on # of files open at any one time.
    //
    for (int i = 0; i < NPROCS; i++)
    {
        if (i == MYPROC)
        {
            Array<int> tranFile = encodeStringForFortran(TransportFile);
            int        lentran  = TransportFile.size();

            FORT_SETUPCKIO(tranFile.dataPtr(), &lentran);
            FORT_INTERPCHEM();
            FORT_CLOSECKIO();
        }
        ParallelDescriptor::Barrier();
    }

    getSpeciesNames();
    getElementNames();
    FORT_GETCKDIMPARAMS(&mMaxreac, &mMaxspec, &mMaxelts,  &mMaxord,
                        &mMaxthrdb, &mMaxtp,  &mMaxsp,    &mMaxspnml);
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
	int len = FORT_GETCKSPECNAME(&i, coded);
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
	int len = FORT_GETCKELTNAME(&i, coded);
	mElementNames[i] = decodeStringFromFortran(coded,len);
    }
    delete [] coded;
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
    FORT_CKINU(&Nids,KI.dataPtr(),&mMaxsp,NU.dataPtr(),&mMaxsp,&fortReacIdx);
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

void
ChemDriver::reactionRateC(FArrayBox&       Cdot,
                             const FArrayBox& C,
                             const FArrayBox& T,
                             Real             Patm,
                             const Box&       box,
                             int              sCompC,
                             int              sCompT,
                             int              sCompCdot) const
{
    const int Nspec = numSpecies();
    BL_ASSERT(Cdot.nComp() >= sCompCdot + Nspec);
    BL_ASSERT(C.nComp() >= sCompC + Nspec);
    BL_ASSERT(T.nComp() > sCompT);

    const Box& mabx = C.box();
    const Box& mbbx = T.box();
    const Box& mobx = Cdot.box();
    
    Box ovlp = box & mabx & mbbx & mobx;
    if( ! ovlp.ok() ) return;
    
    FORT_RRATEC(ovlp.loVect(), ovlp.hiVect(),
                C.dataPtr(sCompC),       ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
                T.dataPtr(sCompT),       ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
                Cdot.dataPtr(sCompCdot), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()),
                &Patm);
}

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
                              const Chem_Evolve solver) const
{
    BL_ASSERT(sCompY+numSpecies() <= Ynew.nComp());
    BL_ASSERT(sCompY+numSpecies() <= Yold.nComp());
    BL_ASSERT(sCompT < Tnew.nComp());
    BL_ASSERT(sCompT < Told.nComp());
    
    BL_ASSERT(Ynew.box().contains(box) && Yold.box().contains(box));
    BL_ASSERT(Tnew.box().contains(box) && Told.box().contains(box));

    if (solver == CKD_Vode)
    {
        FORT_CONPSOLV(box.loVect(), box.hiVect(),
                      Ynew.dataPtr(sCompY), ARLIM(Ynew.loVect()), ARLIM(Ynew.hiVect()),
                      Tnew.dataPtr(sCompT), ARLIM(Tnew.loVect()), ARLIM(Tnew.hiVect()),
                      Yold.dataPtr(sCompY), ARLIM(Yold.loVect()), ARLIM(Yold.hiVect()),
                      Told.dataPtr(sCompT), ARLIM(Told.loVect()), ARLIM(Told.hiVect()),
                      FuncCount.dataPtr(),
                      ARLIM(FuncCount.loVect()), ARLIM(FuncCount.hiVect()),
                      &Patm, &dt);
    }
    else
    {
        FORT_CHEMEQ(box.loVect(), box.hiVect(),
                    Ynew.dataPtr(sCompY), ARLIM(Ynew.loVect()), ARLIM(Ynew.hiVect()),
                    Tnew.dataPtr(sCompT), ARLIM(Tnew.loVect()), ARLIM(Tnew.hiVect()),
                    Yold.dataPtr(sCompY), ARLIM(Yold.loVect()), ARLIM(Yold.hiVect()),
                    Told.dataPtr(sCompT), ARLIM(Told.loVect()), ARLIM(Told.hiVect()),
                    FuncCount.dataPtr(),
                    ARLIM(FuncCount.loVect()), ARLIM(FuncCount.hiVect()),
                    &Patm, &dt);
    }
}

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
    
    FORT_MIXAVG_RHODIFF(box.loVect(), box.hiVect(),
                        rhoD.dataPtr(sCompRD),ARLIM(rhoD.loVect()),ARLIM(rhoD.hiVect()),
                        T.dataPtr(sCompT),    ARLIM(T.loVect()),   ARLIM(T.hiVect()),
                        Y.dataPtr(sCompY),    ARLIM(Y.loVect()),   ARLIM(Y.hiVect()),
                        &Patm);


}

void
ChemDriver::getMixAveragedRhoDiffP(FArrayBox&       rhoD,
                                      const FArrayBox& Y,
                                      const FArrayBox& T,
                                      const FArrayBox& Pmks,
                                      const Box&       box,
                                      int              sCompY,
                                      int              sCompT,
                                      int              sCompP,
                                      int              sCompRD) const
{
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    BL_ASSERT(rhoD.nComp() >= sCompRD+numSpecies());
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Pmks.nComp() > sCompP);
    BL_ASSERT(Y.box().contains(box));
    BL_ASSERT(Pmks.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(rhoD.box().contains(box));
    
    FORT_MIXAVG_RHODIFF_P(box.loVect(), box.hiVect(),
                          rhoD.dataPtr(sCompRD),
                          ARLIM(rhoD.loVect()),ARLIM(rhoD.hiVect()),
                          T.dataPtr(sCompT),   ARLIM(T.loVect()),   ARLIM(T.hiVect()),
                          Y.dataPtr(sCompY),   ARLIM(Y.loVect()),   ARLIM(Y.hiVect()),
                          Pmks.dataPtr(sCompY),
                          ARLIM(Pmks.loVect()),ARLIM(Pmks.hiVect()));
}

void
ChemDriver::getMixCond(FArrayBox&       lambda,
                          const FArrayBox& Y,
                          const FArrayBox& T,
                          const Box&       box,
                          int              sCompY,
                          int              sCompT,
                          int              sCompL) const
{
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    BL_ASSERT(lambda.nComp() >= sCompL);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(lambda.box().contains(box));
    
    FORT_MIX_COND(box.loVect(), box.hiVect(),
                  lambda.dataPtr(sCompL),ARLIM(lambda.loVect()),ARLIM(lambda.hiVect()),
                  T.dataPtr(sCompT),     ARLIM(T.loVect()),     ARLIM(T.hiVect()),
                  Y.dataPtr(sCompY),     ARLIM(Y.loVect()),     ARLIM(Y.hiVect()));
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
ChemDriver::getMixBulkVisc(FArrayBox&       kappa,
                              const FArrayBox& Y,
                              const FArrayBox& T,
                              const Box&       box,
                              int              sCompY,
                              int              sCompT,
                              int              sCompK) const
{
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    BL_ASSERT(kappa.nComp() >= sCompK);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(kappa.box().contains(box));
    
    FORT_MIX_BULK_VISC(box.loVect(), box.hiVect(),
                       kappa.dataPtr(sCompK),ARLIM(kappa.loVect()),ARLIM(kappa.hiVect()),
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
			   int              sCompT) const
{
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(H.nComp() > sCompH);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(H.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    return FORT_TfromHY(box.loVect(), box.hiVect(),
			T.dataPtr(sCompT), ARLIM(T.loVect()), ARLIM(T.hiVect()),
			H.dataPtr(sCompH), ARLIM(H.loVect()), ARLIM(H.hiVect()),
			Y.dataPtr(sCompY), ARLIM(Y.loVect()), ARLIM(Y.hiVect()),
			&mHtoTerrMAX,&mHtoTiterMAX,mTmpData.dataPtr());
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
