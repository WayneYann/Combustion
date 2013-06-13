! ::: 
! ::: ----------------------------------------------------------------
! ::: 
! Passing data from f90 back to C++
subroutine get_method_params(ngrow_c)
  use meth_params_module
  implicit none 
  integer, intent(out) :: ngrow_c
  ngrow_c = NGROW
end subroutine get_method_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 
subroutine set_method_params(dm,Density,Xmom,Eden,Temp,FirstSpec, &
     NUM_STATE, NumSpec, small_dens_in, small_temp_in, small_pres_in, gamma_in)

  use meth_params_module
  use eos_module

  implicit none 

  integer, intent(in) :: dm
  integer, intent(in) :: Density, Xmom, Eden, Temp, FirstSpec, NUM_STATE, NumSpec
  double precision, intent(in) :: small_dens_in, small_temp_in, small_pres_in, gamma_in
  
  integer QLAST

  NVAR = NUM_STATE
  NSPEC = NumSpec

  ! We use these to index into the state "U"
  URHO  = Density   + 1
  UMX   = Xmom      + 1
  if (dm .ge. 2) UMY = UMX + 1
  if (dm .eq. 3) UMZ = UMY + 1
  UEDEN = Eden      + 1
  UTEMP = Temp      + 1
  UFS   = FirstSpec + 1

  ! We use these to index into the state "Q"
  QRHO  = 1
  QU    = 2
  QLAST = 2

  if (dm .ge. 2) then
     QV    = 3
     QLAST = 3
  end if
  
  if (dm .eq. 3) then
     QW    = 4
     QLAST = 4
  end if

  QPRES   = QLAST + 1
  QREINT  = QLAST + 2
  QTEMP   = QLAST + 3

  if (NumSpec .gt. 0) then
     QFS     = QLAST + 4
  else
     QFS     = 0
  end if

  QVAR = QTEMP + NumSpec

  if (small_pres_in > 0.d0) then
     small_pres = small_pres_in
  else
     small_pres = 1.d-50
  end if

  call eos_init(small_dens=small_dens_in, small_temp=small_temp_in, gamma_in=gamma_in)

  call eos_get_small_dens(small_dens)
  call eos_get_small_temp(small_temp)

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
