#ifndef _ChemDriver_H_
#define _ChemDriver_H_

#include <string>

#include <Array.H>

class ChemDriver
{
public:

    ChemDriver (int use_vode);

    ~ChemDriver ();

    void reset();
    
    //
    // Species info.
    //
    int numSpecies () const;
    int numElements () const;
    const Array<std::string>& speciesNames () const;
    const Array<std::string>& elementNames () const;
    int index(const std::string speciesName) const;

    bool isNull;
    bool use_vode;

private:

    void initOnce ();
    
    Array<std::string> mSpeciesNames;
    Array<std::string> mElementNames;
    
    std::string decodeStringFromFortran(const int* coded, int length);
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

#endif
