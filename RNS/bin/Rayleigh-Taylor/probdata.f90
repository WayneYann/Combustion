module probdata_module

  integer, save :: prob_type

  double precision, save :: frac

  double precision, save :: center(3)

  ! RT parameters
  double precision, save :: rho_1, rho_2
  double precision, save :: p0_base
  double precision, save :: L_x

  ! These determine the refinement criteria
  double precision, save :: dengrad
  integer         , save :: max_dengrad_lev

end module probdata_module
