subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use RNS_boundary_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ pamb, v_cf, T_cf, bl_h, max_denerr_lev

  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)
  
  if (namlen .gt. maxlen) then
     write(6,*) 'probin file name too long'
     stop
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! defaults
  pamb = 1.01325d6  ! 1 Patm
  v_cf = 5.5d3
  T_cf = 750.d0
  bl_h = 1.d0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  Twall = T_cf

end subroutine PROBINIT


subroutine rns_initdata(level,time,lo,hi,nscal, &
     state,state_l1,state_l2,state_h1,state_h2, &
     delta,xlo,xhi)

  use probdata_module
  use meth_params_module
  use chemistry_module

  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)
  
  ! local variables
  integer :: i, j, iwrk
  integer :: iN2, iO2
  double precision :: y, rwrk, X0(nspec), Y0(nspec), rho0, e0, T0, state0(NVAR), v

  iN2 = get_species_index("N2")
  iO2 = get_species_index("O2")

  X0 = 0.d0
  X0(iN2) = 0.79d0
  X0(iO2) = 0.21d0

  T0 = T_cf

  call CKXTY(X0, iwrk, rwrk, Y0)
  call CKRHOX(pamb, T0, X0, iwrk, rwrk, rho0)
  call CKUBMS(T0, Y0, iwrk, rwrk, e0)

  state0(URHO)  = rho0
  state0(UMX)   = rho0*v_cf
  state0(UMY)   = 0.d0
  state0(UEDEN) = rho0*(e0+0.5d0*v_cf**2)
  state0(UTEMP) = T0
  state0(UFS:UFS+nspec-1) = rho0*Y0

  !$omp parallel do private(i,j,y,v)
  do j = lo(2), hi(2)
     y = xlo(2) + (j-lo(2)+0.5d0)*delta(2)
     if (y .le. bl_h) then
        v = v_cf * (1.d0-(1.d0-y/bl_h)**50)
     else
        v = v_cf
     end if
     do i = lo(1), hi(1)
        state(i,j,:) = state0
        state(i,j,UMX  ) = rho0*v
        state(i,j,UEDEN) = rho0*(e0 + 0.5d0*v**2)
     end do
  end do
  !$omp end parallel do

end subroutine rns_initdata
