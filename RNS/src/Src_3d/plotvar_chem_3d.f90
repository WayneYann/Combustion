!xxxxxxxxxxxxxxxx todo

subroutine rns_makeplotvar(lo, hi, dx, &
     prim, pm_l1, pm_l2, pm_l3, pm_h1, pm_h2, pm_h3, &
     plot, pt_l1, pt_l2, pt_l3, pt_h1, pt_h2, pt_h3, &
     npv, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, &
     icomp_X, icomp_omegadot, icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, &
     fuelID)
  use meth_params_module, only : NVAR, QRHO, QU, QV, QW, QTEMP, QFY, NSPEC
  use eos_module, only : eos_get_c
  implicit none

  double precision, intent(in) :: dx(3)
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: pm_l1, pm_h1, pm_l2, pm_h2, pm_l3, pm_h3
  integer, intent(in) :: pt_l1, pt_h1, pt_l2, pt_h2, pt_l3, pt_h3
  integer, intent(in) :: npv, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, &
       icomp_X, icomp_omegadot,	icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, &
       fuelID
  double precision,intent(in   )::prim(pm_l1:pm_h1,pm_l2:pm_h2,pm_l3:pm_h3,NVAR) 
  double precision,intent(inout)::plot(pt_l1:pt_h1,pt_l2:pt_h2,pt_l3:pt_h3,0:npv-1) 
    
  integer :: i, j, k
  double precision :: Yt(NSPEC), v, c, dxinv(3)


end subroutine rns_makeplotvar
