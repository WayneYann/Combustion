#include "ChemDriver.H"
#include "ChemDriver_F.H"

#include <ParmParse.H>

namespace
{
    bool initialized = false;

    void ChemDriver_Finalize() { initialized = false; }
}

void ChemDriver::reset()
{
    initialized = false;
}

ChemDriver::ChemDriver (int use_vode_in)
{
    isNull = false;
    use_vode = use_vode_in;

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
    int max_len_elem_name, max_len_spec_name, nelem, nspec;
    
    BL_FORT_PROC_CALL(CD_INITCHEM, cd_initchem) 
	(max_len_elem_name, max_len_spec_name, nelem, nspec);
    
    mElementNames.clear();
    mElementNames.resize(nelem);
    
    mSpeciesNames.clear();
    mSpeciesNames.resize(nspec);
    
    int* elem_name_f = new int[max_len_elem_name];
    int* spec_name_f = new int[max_len_spec_name];
    
    for (int i=0; i<nelem; i++) 
    {
	int ielem = i+1;
	int lname;
	BL_FORT_PROC_CALL(CD_GETELEMNAME, cd_getelemname)
	    (ielem, lname, elem_name_f);
	mElementNames[i] = decodeStringFromFortran(elem_name_f, lname);
    }
    
    for (int i=0; i<nspec; i++) 
    {
	int ispec = i+1;
	int lname;
	BL_FORT_PROC_CALL(CD_GETSPECNAME, cd_getspecname)
	    (ispec, lname, spec_name_f);
	mSpeciesNames[i] = decodeStringFromFortran(spec_name_f, lname);
    }
    
    delete [] elem_name_f;
    delete [] spec_name_f;
    
    // vode
    if (use_vode)
    {
	int  itol = 1;
	Real rtol = 1.e-10;
	Real atol = 1.e-10;
	int  order = 100;  // vode will reduce it to a lower value
	int  maxstep = 2000;  // max step vode takes
	int  use_ajac = 1; // use analytic Jacobian
	int  save_ajac = 1; // reuse analytic Jacobian
	int  always_new_j = 0;
	int  stiff = 1;
	int  verbose = 0;
	
	ParmParse ppv("vode");
	ppv.query("itol", itol);
	ppv.query("rtol", rtol);
	ppv.query("atol", atol);
	ppv.query("order", order);
	ppv.query("maxstep", maxstep);
	ppv.query("use_ajac", use_ajac);
	ppv.query("save_ajac", save_ajac);
	ppv.query("always_new_j", always_new_j);
	ppv.query("stiff", stiff);
	ppv.query("v", verbose);
	ppv.query("verbose", verbose);
	
	BL_ASSERT(rtol > 0);
	BL_ASSERT(atol > 0);
	BL_ASSERT(itol == 1 || itol == 2);
    
	int neq = nspec+1; 
	BL_FORT_PROC_CALL(CD_INITVODE, cd_initvode)
	    (neq, verbose, itol, rtol, atol, order, 
	     maxstep, use_ajac, save_ajac, always_new_j, stiff); 
    }
    else
    {
	Real rtol = 1.e-10;
	Real atol = 1.e-10;
	int  order = 6;  
	int  verbose = 0;
	int reuse_jac = 1;
	int multipoint = 1;

	ParmParse ppb("bdf");
	ppb.query("rtol", rtol);
	ppb.query("atol", atol);
	ppb.query("order", order);
	ppb.query("v", verbose);
	ppb.query("verbose", verbose);
	ppb.query("reuse_jac", reuse_jac); 
	ppb.query("multipoint", multipoint);

	int npt = 1;
	if (multipoint)
	{
#if (BL_SPACEDIM == 1)
	    npt = 2;
#elif (BL_SPACEDIM == 2)
	    npt = 4;
#else
	    npt = 8;
#endif
	}

	int neq = nspec+1; 
	BL_FORT_PROC_CALL(CD_INITBDF, cd_initbdf)
	    (neq, npt, verbose, rtol, atol, order, reuse_jac);
    }


    // eglib
    int use_bulk_visc = 1;

    ParmParse ppe("eglib");
    ppe.query("use_bulk_visc", use_bulk_visc);
    
    BL_FORT_PROC_CALL(CD_INITEGLIB, cd_initeglib)
	(use_bulk_visc);
}


ChemDriver::~ChemDriver ()
{
    BL_FORT_PROC_CALL(CD_CLOSECHEM, cd_closechem)();
    if (use_vode)
    {
	BL_FORT_PROC_CALL(CD_CLOSEVODE, cd_closevode)();
    }
    else
    {
	BL_FORT_PROC_CALL(CD_CLOSEBDF, cd_closebdf)();
    }
    BL_FORT_PROC_CALL(CD_CLOSEEGLIB, cd_closeeglib)();
}


std::string
ChemDriver::decodeStringFromFortran(const int* coded, int length)
{
    std::string result;
    for (int i = 0; i < length; ++i)
        result += coded[i];
    return result;
}


int
ChemDriver::index(const std::string speciesName) const
{
    for (int i = 0; i < mSpeciesNames.size(); i++)
    {
	if (speciesName == mSpeciesNames[i])
	{
	    return i;
	}
    }
    return -1;
}

