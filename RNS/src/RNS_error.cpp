#include <winstd.H>

#include "RNS.H"
#include "RNS_F.H"

using std::string;

typedef StateDescriptor::BndryFunc BndryFunc;

void
RNS::ErrorSetUp ()
{
    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //

    // This routine uses the special evaluation of the second derivative
    //   and can be called with any variable.  The quantity "lap_var"
    //   must first be defined as a derived quantity.  In this example
    //   "lap_var" is the weighted Laplacian of density.
//  derive_lst.add("lap_var",IndexType::TheCellType(),1,FORT_DERLAPVAR,grow_box_by_two);
//  derive_lst.addComponent("lap_var",desc_lst,State_Type,Density,1);
//  err_list.add("lap_var",1,ErrorRec::Special,
//		 BL_FORT_PROC_CALL(RNS_SPECIAL_ERROR,cns_special_error));

    err_list.add("density",1,ErrorRec::Special,
		 BL_FORT_PROC_CALL(RNS_DENERROR,cns_denerror));
//  err_list.add("Temp",1,ErrorRec::Special,
//		 BL_FORT_PROC_CALL(RNS_TEMPERROR,cns_temperror));
    err_list.add("pressure",1,ErrorRec::Special,
		 BL_FORT_PROC_CALL(RNS_PRESSERROR,cns_presserror));
    err_list.add("x_velocity",1,ErrorRec::Special,
                   BL_FORT_PROC_CALL(RNS_VELERROR,cns_velerror));
#if (BL_SPACEDIM >= 2)
    err_list.add("y_velocity",1,ErrorRec::Special,
                   BL_FORT_PROC_CALL(RNS_VELERROR,cns_velerror));
#endif
#if (BL_SPACEDIM == 3)
    err_list.add("z_velocity",1,ErrorRec::Special,
                   BL_FORT_PROC_CALL(RNS_VELERROR,cns_velerror));
#endif

//   err_list.add("entropy",1,ErrorRec::Special,
//		 BL_FORT_PROC_CALL(RNS_ENTERROR,cns_enterror));
}
