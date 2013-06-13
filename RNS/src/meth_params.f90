module meth_params_module

  implicit none

  integer, parameter     :: NGROW = 3  ! fifth-order WENO

  ! conserved variables
  integer         , save :: NVAR, NSPEC
  integer         , save :: URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS

  ! primitive variables
  integer         , save :: QVAR
  integer         , save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QFS

  double precision, save :: small_dens, small_temp, small_pres  

end module meth_params_module
