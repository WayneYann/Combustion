#ifndef _ChemDriver_H_
#define _ChemDriver_H_

#include <string>

#include <Array.H>
#include <FArrayBox.H>

class ChemDriver
{
public:

    ChemDriver (const std::string& TransportFile);
		  
    ~ChemDriver () {}
    
    enum Chem_Evolve { CKD_Vode, CKD_ChemEQ2 };

    void solveTransient (FArrayBox&        Ynew,
                         FArrayBox&        Tnew,
                         const FArrayBox&  Yold,
                         const FArrayBox&  Told,
                         FArrayBox&        FuncCount,
                         const Box&        box,
                         int               sCompY,
                         int               sCompT,
                         Real              dt,
                         Real              Patm,
                         const Chem_Evolve solver=CKD_Vode,
                         FArrayBox*        chemDiag=0) const;

    void set_verbose_vode ();
    void set_max_vode_subcycles (int max_cyc);
    void set_species_Yscales (const std::string& scalesFile);
    //
    // Species info.
    //
    int numSpecies () const;
    int numElements () const;
    int numReactions () const;
    const Array<std::string>& speciesNames () const;
    const Array<std::string>& elementNames () const;
    Array<Real> speciesMolecWt () const;
    Array<Real> elementAtomicWt () const;
    int index (const std::string& specName) const;
    int indexElt (const std::string& eltName) const;
    Array<int> reactionsWithXonL (const std::string& lSpecName) const;
    Array<int> reactionsWithXonR (const std::string& rSpecName) const;
    Array<std::pair<std::string,int> > specCoeffsInReactions (int reacIdx) const;
    std::string reactionString (int reacIdx) const;

    int numberOfElementXinSpeciesY (const std::string& eltX,
                                    const std::string& spcY) const;
    //
    // Thermo info.
    //
    void getRhoGivenPTY (FArrayBox&       Rho,
                         Real             Patm,
                         const FArrayBox& T,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompT,
                         int              sCompY,
                         int              sCompR) const;

    void getRhoGivenPvTY (FArrayBox&       Rho,
                          const FArrayBox& P,
                          const FArrayBox& T,
                          const FArrayBox& Y,
                          const Box&       box,
                          int              sCompP,
                          int              sCompT,
                          int              sCompY,
                          int              sCompR) const;

    void getPGivenRTY (FArrayBox&       p,
                       const FArrayBox& Rho,
                       const FArrayBox& T,
                       const FArrayBox& Y,
                       const Box&       box,
                       int              sCompR,
                       int              sCompT,
                       int              sCompY,
                       int              sCompP) const;

    void getTGivenPRY (FArrayBox&       T,
                       Real             p,
                       const FArrayBox& Rho,
                       const FArrayBox& Y,
                       const Box&       box,
                       int              sCompR,
                       int              sCompY,
                       int              sCompT) const;

    void getCpmixGivenTY (FArrayBox&       cpmix,
                          const FArrayBox& T,
                          const FArrayBox& Y,
                          const Box&       box,
                          int              sCompT,
                          int              sCompY,
                          int              sCompCp) const;

    void getCvmixGivenTY (FArrayBox&       cvmix,
                          const FArrayBox& T,
                          const FArrayBox& Y,
                          const Box&       box,
                          int              sCompT,
                          int              sCompY,
                          int              sCompCv) const;

    void getHmixGivenTY (FArrayBox&       hmix,
                         const FArrayBox& T,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompT,
                         int              sCompY,
                         int              sCompH) const;

    void getMwmixGivenY (FArrayBox&       mwmix,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompY,
                         int              sCompMw) const;

    void getCpGivenT (FArrayBox&       cp,
                      const FArrayBox& T,
                      const Box&       box,
                      int              sCompT,
                      int              sCompCp) const;

    void getHGivenT (FArrayBox&       h,
                     const FArrayBox& T,
                     const Box&       box,
                     int              sCompT,
                     int              sCompH) const;

    Real getRuniversal () const;
    Real getP1atm_MKS () const;
    //
    // Compute T that satisfies hmix=sum(h(T)), returns max Newton iterations
    // on any point over grid, return -1 if T jumped out of bounds during solve,
    // and -2 if solve failed anywhere.  Note that a temporary is not used, and
    // the solver kicks out when/if it fails, so it may exit after converting
    // only part of the temperature array.  Save a copy of T and check the return
    // code if you want to be extra careful...
    //
    int getTGivenHY (FArrayBox&       T,
                     const FArrayBox& H,
                     const FArrayBox& Y,
                     const Box&       box,
                     int              sCompH,
                     int              sCompY,
                     int              sCompT) const;

    Array<Real> massFracToMoleFrac (const Array<Real>& Y) const;
    Array<Real> moleFracToMassFrac (const Array<Real>& X) const;

    void molarProduction (FArrayBox&       Q,
                          const std::string&   specName,
                          const FArrayBox& C,
                          const FArrayBox& T,
                          const Box&       box,
                          int              sCompC,
                          int              sCompT,
                          int              sCompQ) const;
    //
    // Compute heat release (J/(m^3.s)) based on the temp, press and mass fractions
    //
    void heatRelease (FArrayBox&       Q,
                      const FArrayBox& Y,
                      const FArrayBox& T,
                      Real             Patm,
                      const Box&       box,
                      int              sCompY,
                      int              sCompT,
                      int              sCompQ) const;
    //
    // Compute dY/dt based on the input temp, press and mass fractions.
    //
    void reactionRateY (FArrayBox&       Ydot,
                        const FArrayBox& Y,
                        const FArrayBox& T,
                        Real             Patm,
                        const Box&       box,
                        int              sCompY,
                        int              sCompT,
                        int              sCompYdot) const;

    void reactionRateC (FArrayBox&       Cdot,
                        const FArrayBox& C,
                        const FArrayBox& T,
                        Real             Patm,
                        const Box&       box,
                        int              sCompC,
                        int              sCompT,
                        int              sCompCdot) const;

    void fwdRevReacRatesGivenXTP (FArrayBox&        FwdK,
                                  FArrayBox&        RevK,
                                  const Array<int>& rxnIDs,
                                  const FArrayBox&  X,
                                  const FArrayBox&  T,
                                  Real              Patm,
                                  const Box&        box,
                                  int               sCompX,
                                  int               sCompT,
                                  int               sCompFwdK,
                                  int               sCompRevK) const;
    
    void massFracToMoleFrac (FArrayBox&       X,
                             const FArrayBox& Y,
                             const Box&       box,
                             int              sCompY,
                             int              sCompX) const;

    void moleFracToMassFrac (FArrayBox&       Y,
                             const FArrayBox& X,
                             const Box&       box,
                             int              sCompX,
                             int              sCompY) const;

    void massFracToMolarConc (FArrayBox&       C,
                              const FArrayBox& Y,
                              const FArrayBox& T,
                              Real             Patm,
                              const Box&       box,
                              int              sCompY,
                              int              sCompT,
                              int              sCompC) const;

    void massFracToMolarConc (FArrayBox&       C,
                              const FArrayBox& Y,
                              const FArrayBox& T,
                              const FArrayBox& Rho,
                              const Box&       box,
                              int              sCompR,
                              int              sCompY,
                              int              sCompT,
                              int              sCompC) const;

    void molarConcToMoleFrac (FArrayBox&       X,
                              const FArrayBox& C,
                              const Box&       box,
                              int              sCompC,
                              int              sCompX) const;
    //
    // Normalize mass fractions to prevent negatives.
    //
    void normalizeMassFrac (FArrayBox&       Ynorm,
                            const FArrayBox& Y,
                            const std::string&   excessSpecies,
                            const Box&       box,
                            int              sCompY,
                            int              sCompYnorm) const;
    //
    // Chemical Diffusivities.
    //
    void getMixAveragedRhoDiff (FArrayBox&       rhoD,
                                const FArrayBox& Y,
                                const FArrayBox& T,
                                Real             Patm,
                                const Box&       box,
                                int              sCompY,
                                int              sCompT,
                                int              sCompRD) const;
    
    void getMixAveragedRhoDiffP (FArrayBox&       rhoD,
                                 const FArrayBox& Y,
                                 const FArrayBox& T,
                                 const FArrayBox& Pmks,
                                 const Box&       box,
                                 int              sCompY,
                                 int              sCompT,
                                 int              sCompP,
                                 int              sCompRD) const;
    //
    // Thermal conductivity.
    //
    void getMixCond (FArrayBox&       lambda,
                     const FArrayBox& Y,
                     const FArrayBox& T,
                     const Box&       box,
                     int              sCompY,
                     int              sCompT,
                     int              sCompL) const;
    //
    // Viscosity.
    //
    void getMixShearVisc (FArrayBox&       eta,
                          const FArrayBox& Y,
                          const FArrayBox& T,
                          const Box&       box,
                          int              sCompY,
                          int              sCompT,
                          int              sCompE) const;
    
    void getMixBulkVisc (FArrayBox&       kappa,
                         const FArrayBox& Y,
                         const FArrayBox& T,
                         const Box&       box,
                         int              sCompY,
                         int              sCompT,
                         int              sCompK) const;
    
    void getElementMoles (FArrayBox&       C_elt,
                          const std::string&   name,
                          const FArrayBox& C,
                          const Box&       box,
                          int              sCompC,
                          int              sCompC_elt) const;
    //
    // Optically thin radiation model.
    //
    void getOTradLoss_TDF (FArrayBox&       Qloss,
                           const FArrayBox& T,
                           const FArrayBox& X,
                           const Real       Patm,
                           const Real       T_bg,
                           const Box&       box,
                           int              sCompX,
                           int              sCompT,
                           int              sCompQ) const;
    //
    // H - to - T solve parameter access.
    //
    Real getHtoTerrMAX () const;
    void setHtoTerrMAX (Real err);
    int getHtoTiterMAX () const;
    void setHtoTiterMAX (int err);    
    //
    // Handy functions.
    //
    static Array<int> encodeStringForFortran(const std::string& astr);
    static std::string decodeStringFromFortran(const int* coded, int length);

protected:

    void getSpeciesNames ();
    void getElementNames ();

private:

    void initOnce (const std::string& TransportFile);
    
    static bool mInitialized;

    Array<std::string> mSpeciesNames;
    Array<std::string> mElementNames;
    Real               mHtoTerrMAX;
    int                mHtoTiterMAX;
    Array<Real>        mTmpData;
    int mMaxreac, mMaxspec, mMaxelts, mMaxord, mMaxthrdb, mMaxtp, mMaxsp, mMaxspnml;
};

inline
int
ChemDriver::numSpecies () const
{
    return mSpeciesNames.size();
}

inline
int
ChemDriver::numElements () const
{
    return mElementNames.size();
}

inline
const Array<std::string>&
ChemDriver::speciesNames () const
{
    return mSpeciesNames;
}

inline
const Array<std::string>&
ChemDriver::elementNames () const
{
    return mElementNames;
}

inline
int
ChemDriver::index (const std::string& specName) const
{
    for (int i=0; i<mSpeciesNames.size(); i++)
    {
        if (specName == mSpeciesNames[i])
            return i;
    }
    return -1;
}

inline
int
ChemDriver::indexElt (const std::string& eltName) const
{
    for (int i=0; i<mElementNames.size(); i++)
    {
        if (eltName == mElementNames[i])
            return i;
    }
    return -1;
}

inline
Real
ChemDriver::getHtoTerrMAX () const
{
    return mHtoTerrMAX;
}

inline
void
ChemDriver::setHtoTerrMAX (Real err)
{
    mHtoTerrMAX = err;
}

inline
int
ChemDriver::getHtoTiterMAX () const
{
    return mHtoTiterMAX;
}

inline
void
ChemDriver::setHtoTiterMAX (int iter)
{
    if (iter != mHtoTiterMAX)
    {
        mHtoTiterMAX = iter;
        mTmpData.resize(mHtoTiterMAX);
    }
}

#endif /*_ChemDriver_H_*/
