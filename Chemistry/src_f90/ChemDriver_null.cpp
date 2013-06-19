// Null Chemistry

#include "ChemDriver.H"

namespace
{
    bool initialized = false;
    
    void ChemDriver_Finalize() { initialized = false; }
}

ChemDriver::ChemDriver ()
{
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
    mSpeciesNames.resize(1);
    mSpeciesNames[0] = "X";

    mElementNames.clear();
    mElementNames.resize(1);
    mElementNames[0] = "X";
}


ChemDriver::~ChemDriver ()
{
    ;
}
