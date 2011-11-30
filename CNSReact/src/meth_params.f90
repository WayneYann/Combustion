
! This module stores the runtime parameters.  
! These parameter are initialized in set_method_params().

module meth_params_module

  implicit none

  double precision, save :: difmag        ! used only in consup to weight the divu contributin
  integer         , save :: iorder        ! used only in uslope and uflaten

  integer, parameter     :: NHYP    = 4
  integer, parameter     :: MAXADV  = 2

  ! NTHERM: number of thermodynamic variables
  integer         , save :: NTHERM, NVAR
  integer         , save :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFA, UFS

  ! QTHERM: number of primitive variables
  integer         , save :: QTHERM, QVAR
  integer         , save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP
  integer         , save :: QFA, QFS

  integer         , save :: nadv

  double precision, save :: small_dens, small_temp, small_pres  

  integer         , save :: allow_negative_energy
  integer         , save :: ppm_type
  integer         , save :: normalize_species

end module meth_params_module
