module plotvar_index_module

  implicit none

  ! the total number of plot components
  integer, save :: n_plot_comps = 0
  integer, save :: icomp_rho=0, icomp_vel=0, icomp_pres=0, icomp_temp=0, icomp_eint=0, &
       icomp_h=0, icomp_rhoh=0, icomp_cs, icomp_magvel, icomp_Mach, &
       icomp_divu=0, icomp_magvort=0, icomp_Y=0, icomp_X=0, &
       icomp_wbar=0, icomp_hspec=0, icomp_omegadot=0, icomp_dYdt=0, icomp_heatRelease=0, &
       icomp_fuelConsumption=0
  ! "burn" includes omegadot, dYdt, heatRelease & fuelconsumption
  integer, save :: icomp_burn=0 
  integer, save :: ifuel = -1
  integer, save :: nburn = 0
  integer, save :: ib_omegadot=0, ib_dYdt=0, ib_heatRelease=0, &
       ib_fuelConsumption=0

end module plotvar_index_module

