#include "ChemDriver.H"
#include "ChemDriver_F.H"

#include <iostream>
using namespace std;

ChemDriver::ChemDriver ()
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
}


ChemDriver::~ChemDriver ()
{
  BL_FORT_PROC_CALL(CD_closeCHEM, cd_closechem)();
}


std::string
ChemDriver::decodeStringFromFortran(const int* coded, int length)
{
    std::string result;
    for (int i = 0; i < length; ++i)
        result += coded[i];
    return result;
}
