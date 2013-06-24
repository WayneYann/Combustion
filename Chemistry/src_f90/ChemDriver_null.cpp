// Null Chemistry

#include "ChemDriver.H"

namespace
{
    bool initialized = false;
    
    void ChemDriver_Finalize() { initialized = false; }
}

ChemDriver::ChemDriver ()
{
    isNull = true;

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
