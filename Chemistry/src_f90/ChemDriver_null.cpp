// Null Chemistry

#include "ChemDriver.H"

namespace
{
    bool initialized = false;
    
    void ChemDriver_Finalize() { initialized = false; }
}

ChemDriver::ChemDriver (int use_vode_in)
{
    isNull = true;
    use_vode = false;

    if (!initialized) 
    {
	initOnce();
	BoxLib::ExecOnFinalize(ChemDriver_Finalize);
	initialized = true;
    }
}

void 
ChemDriver::initOnce ()
{
    mSpeciesNames.clear();
    mSpeciesNames.resize(2);
    mSpeciesNames[0] = "X";
    mSpeciesNames[1] = "Y";

    mElementNames.clear();
    mElementNames.resize(2);
    mElementNames[0] = "X";
    mElementNames[1] = "Y";
}


ChemDriver::~ChemDriver ()
{
    ;
}

int
ChemDriver::index(const std::string speciesName) const
{
    return -1;
}
