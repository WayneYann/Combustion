// Null Chemistry

#include "ChemDriver.H"

namespace
{
    bool initialized = false;

    void ChemDriver_Finalize() { initialized = false; }
}

ChemDriver::ChemDriver ()
{
  if (!initialized) {
    initOnce();
    BoxLib::ExecOnFinalize(ChemDriver_Finalize);
    initialized = true;
  }
}

void 
ChemDriver::initOnce ()
{
  ;
}


ChemDriver::~ChemDriver ()
{
  ;
}
