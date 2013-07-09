
#include "RNS.H"
#include "RNS_F.H"

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
//		 BL_FORT_PROC_CALL(RNS_SPECIAL_ERROR,rns_special_error));

    if (do_density_ref)
    {
	err_list.add("density",1,ErrorRec::Special,
		     BL_FORT_PROC_CALL(RNS_DENERROR,rns_denerror));
    }

    if (do_temperature_ref)
    {
	err_list.add("Temp",1,ErrorRec::Special,
		     BL_FORT_PROC_CALL(RNS_TEMPERROR,rns_temperror));
    }

    if (do_pressure_ref)
    {
	err_list.add("pressure",1,ErrorRec::Special,
		     BL_FORT_PROC_CALL(RNS_PRESSERROR,rns_presserror));
    }

    if (do_velocity_ref)
    {
	err_list.add("x_velocity",1,ErrorRec::Special,
		     BL_FORT_PROC_CALL(RNS_VELERROR,rns_velerror));
#if (BL_SPACEDIM >= 2)
	err_list.add("y_velocity",1,ErrorRec::Special,
		     BL_FORT_PROC_CALL(RNS_VELERROR,rns_velerror));
#endif
#if (BL_SPACEDIM == 3)
	err_list.add("z_velocity",1,ErrorRec::Special,
		     BL_FORT_PROC_CALL(RNS_VELERROR,rns_velerror));
#endif
    }

    if (do_vorticity_ref)
    {
#if (BL_SPACEDIM >= 2)
	err_list.add("magvort",1,ErrorRec::Special,
		     BL_FORT_PROC_CALL(RNS_VORTERROR,rns_vorterror));
#endif
    }

}
