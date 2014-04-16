! ::: 
! ::: ----------------------------------------------------------------
! ::: 
! Passing data from f90 back to C++
subroutine get_method_params(ngrow_c, nriemann_c)
  use meth_params_module
  implicit none 
  integer, intent(out) :: ngrow_c, nriemann_c
  ngrow_c = NGROW
  nriemann_c = nriemann
end subroutine get_method_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 
subroutine set_method_params(dm,Density,Xmom,Eden,Temp,FirstSpec, &
     NUM_STATE, NumSpec, small_dens_in, small_temp_in, small_pres_in, &
     gamma_in, grav_in, Tref_in, riemann_in, difmag_in, blocksize, &
     do_weno_in, do_quadrature_weno_in, do_comp_weno_in, &
     use_vode_in, do_cc_burning_in, split_burning_in, do_BE_burning_in, new_J_cell_in)

  use meth_params_module
  use eos_module

  implicit none 

  integer, intent(in) :: dm
  integer, intent(in) :: Density, Xmom, Eden, Temp, FirstSpec, NUM_STATE, NumSpec, &
       riemann_in, blocksize(*), do_weno_in, do_quadrature_weno_in, do_comp_weno_in, &
       use_vode_in, do_cc_burning_in, split_burning_in, do_BE_burning_in, new_J_cell_in
  double precision, intent(in) :: small_dens_in, small_temp_in, small_pres_in, &
       gamma_in, grav_in, Tref_in, difmag_in
  
  ndim = dm

!  NCHARV = dm + 1 + NumSpec ! momentum + energy + rhoY
!  CFS = dm + 2
! We always do characteristic decomposition in 3D
  NCHARV = 3 + 1 + NumSpec ! momentum + energy + rhoY
  CFS = 3 + 2

  NVAR = NUM_STATE
  NSPEC = NumSpec

  ! We use these to index into the state "U"
  URHO  = Density   + 1
  UMX   = Xmom      + 1
  if (dm .ge. 2) then
     UMY = UMX + 1
  else
     UMY = 0
  end if
  if (dm .eq. 3) then
     UMZ = UMY + 1
  else
     UMZ = 0
  end if
  UEDEN = Eden      + 1
  UTEMP = Temp      + 1
  UFS   = FirstSpec + 1

  ! We use these to index into the state "Q"
  QRHO  = URHO
  QU    = UMX
  QV    = UMY
  QW    = UMZ
  QPRES = UEDEN
  QTEMP = UTEMP
  QFY   = UFS
  QFX   = UFS +   NumSpec
  QFH   = UFS + 2*NumSpec

  QCVAR = QTEMP + 2*NumSpec
  QFVAR = QTEMP + 3*NumSpec

  if (small_pres_in > 0.d0) then
     small_pres = small_pres_in
  else
     small_pres = 1.d-50
  end if

  call eos_init(small_dens=small_dens_in, small_temp=small_temp_in, &
       gamma_in=gamma_in, Tref_in=Tref_in)

  call eos_get_small_dens(small_dens)
  call eos_get_small_temp(small_temp)

  gravity = grav_in

  riemann_solver = riemann_in
  difmag = difmag_in

  xblksize = blocksize(1)
  if (dm .ge. 2) then
     yblksize = blocksize(2)
  end if
  if (dm .ge. 3) then
     zblksize = blocksize(3)
  end if

  nthreads=0
  !$omp parallel reduction(+:nthreads)
  nthreads = nthreads + 1
  !$omp end parallel

  do_weno = (do_weno_in .ne. 0)
  do_quadrature_weno = (do_quadrature_weno_in .ne. 0)
  do_component_weno = (do_comp_weno_in .ne. 0)

  use_vode = (use_vode_in .ne. 0)
  do_cc_burning = (do_cc_burning_in .ne. 0)
  split_burning = (split_burning_in .ne. 0)
  do_BE_burning = (do_BE_burning_in .ne. 0)
  new_J_cell    = (new_J_cell_in    .ne. 0)

end subroutine set_method_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

subroutine set_problem_params(dm,physbc_lo_in,physbc_hi_in, phys_prob_lo_in,   &
     phys_prob_hi_in, Outflow_in,Symmetry_in,coord_type_in)

  use prob_params_module

  implicit none 
  
  integer, intent(in) :: dm
  integer, intent(in) :: physbc_lo_in(dm),physbc_hi_in(dm)
  integer, intent(in) :: Outflow_in
  integer, intent(in) :: Symmetry_in
  integer, intent(in) :: coord_type_in
  double precision, intent(in) :: phys_prob_lo_in(dm),phys_prob_hi_in(dm)

  allocate(physbc_lo(dm))
  allocate(physbc_hi(dm))
  allocate(phys_prob_lo(dm))
  allocate(phys_prob_hi(dm))
  
  physbc_lo(:) = physbc_lo_in(:)
  physbc_hi(:) = physbc_hi_in(:)
  
  phys_prob_lo(:) = phys_prob_lo_in(:)
  phys_prob_hi(:) = phys_prob_hi_in(:)
  
  Outflow  = Outflow_in
  Symmetry = Symmetry_in
  
  coord_type = coord_type_in

end subroutine set_problem_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 
subroutine rns_set_special_tagging_flag(dummy,flag) 
  use probdata_module
  implicit none
  double precision :: dummy 
  integer          :: flag
end subroutine rns_set_special_tagging_flag
