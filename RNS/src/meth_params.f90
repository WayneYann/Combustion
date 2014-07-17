module meth_params_module

  implicit none

  ! fifth-order WENO needs 3 ghost cells
  integer, parameter :: NGROW = 3  

  ! Riemann solvers
  integer, save :: riemann_solver
  integer, parameter :: HLL_solver = 0
  integer, parameter :: JBB_solver = 1
  integer, parameter :: HLLC_solver = 2
  integer, parameter :: nriemann   = 3

  double precision, save :: difmag, HLL_factor
  
  integer, save :: ndim   ! spatial dimension
  integer, save :: NCHARV ! number of characteristic variables
  integer, save :: CFS    ! first species in characteristic variables

  ! conserved variables
  integer         , save :: NVAR, NSPEC
  integer         , save :: URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS

  ! primitive variables
  integer         , save :: QCVAR, QFVAR  ! cell-centered and face
  integer         , save :: QRHO, QU, QV, QW, QPRES, QTEMP, QFY, QFX, QFH

  double precision, save :: small_dens, small_temp, small_pres  

  double precision, save :: gravity

  ! blocking, values will be reset in set_method_params
  integer, save :: xblksize=2048, yblksize=2048, zblksize=2048
  integer, save :: nthreads=1

  logical, save :: do_weno, do_quadrature_weno, do_component_weno;

  ! chemistry
  logical, save :: use_vode
  logical, save :: new_J_cell
  integer, save :: chem_solver
  integer, parameter :: cc_burning = 0
  integer, parameter :: Gauss_burning = 1
  integer, parameter :: split_burning = 2
  integer, parameter :: BEcc_burning = 3
  integer, parameter :: BEGp_burning = 4
  integer, parameter :: nchemsolver = 5
  logical, save :: chem_do_weno

end module meth_params_module
