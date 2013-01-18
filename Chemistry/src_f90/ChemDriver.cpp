#include "ChemDriver.H"
#include "ChemDriver_F.H"

#include <ParmParse.H>

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
  int max_len_elem_name, max_len_spec_name, nelem, nspec;

  BL_FORT_PROC_CALL(CD_INITCHEM, cd_initchem) 
    (max_len_elem_name, max_len_spec_name, nelem, nspec);

  mElementNames.clear();
  mElementNames.resize(nelem);

  mSpeciesNames.clear();
  mSpeciesNames.resize(nspec);

  int* elem_name_f = new int[max_len_elem_name];
  int* spec_name_f = new int[max_len_spec_name];

  for (int i=0; i<nelem; i++) {
    int ielem = i+1;
    int lname;
    BL_FORT_PROC_CALL(CD_GETELEMNAME, cd_getelemname)
      (ielem, lname, elem_name_f);
    mElementNames[i] = decodeStringFromFortran(elem_name_f, lname);
  }

  for (int i=0; i<nspec; i++) {
    int ispec = i+1;
    int lname;
    BL_FORT_PROC_CALL(CD_GETSPECNAME, cd_getspecname)
      (ispec, lname, spec_name_f);
    mSpeciesNames[i] = decodeStringFromFortran(spec_name_f, lname);
  }

  delete [] elem_name_f;
  delete [] spec_name_f;

  // vode
  int  itol = 1;
  Real rtol = 1.e-10;
  Real atol = 1.e-10;

  ParmParse pp("vode");
  pp.query("itol", itol);
  pp.query("rtol", rtol);
  pp.query("atol", atol);

  BL_ASSERT(rtol > 0);
  BL_ASSERT(atol > 0);
  BL_ASSERT(itol == 1 || itol == 2);

  int neq = nspec; 
  BL_FORT_PROC_CALL(CD_INITVODE, cd_initvode)
    (neq, itol, rtol, atol);  
}


ChemDriver::~ChemDriver ()
{
  BL_FORT_PROC_CALL(CD_CLOSECHEM, cd_closechem)();
  BL_FORT_PROC_CALL(CD_CLOSEVODE, cd_closevode)();
}


std::string
ChemDriver::decodeStringFromFortran(const int* coded, int length)
{
    std::string result;
    for (int i = 0; i < length; ++i)
        result += coded[i];
    return result;
}
