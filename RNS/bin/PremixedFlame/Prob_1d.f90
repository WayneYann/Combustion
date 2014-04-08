
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use chemistry_module, only : Patm, nspecies, get_species_index
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(1), probhi(1)

  integer untin,i
  integer :: iH2, iO2, iN2, iH2O, iHO2, iO, iH2O2, iH, iOH

  namelist /fortin/ init_flm_position, init_flm_width

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
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
         
! set namelist defaults

  init_flm_position = 0.d0
  init_flm_width = 0.01d0

!     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  Length(1) = probhi(1) - problo(1)

  

end subroutine PROBINIT

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine rns_initdata(level,time,lo,hi,nscal, &
     state,state_l1,state_h1,delta,xlo,xhi)

  use eos_module, only : eos_given_PTX
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC
  use chemistry_module, only : Patm, nspecies, get_species_index

  implicit none
  integer level, nscal
  integer lo(1), hi(1)
  integer state_l1,state_h1
  double precision state(state_l1:state_h1,NVAR)
  double precision time, delta(1)
  double precision xlo(1), xhi(1)
  
  integer :: i, n
  integer, save :: iH2=-100, iO2, iN2, iH2O, iHO2, iO, iH2O2, iH, iOH
  double precision :: Xt(nspec), Yt(nspec), Yti(nspecies)
  double precision :: xcen, rhot,Tt,et,Pt, eta, eta01, eta00

  if (nspecies .ne. NSPEC) then
     write(6,*)"nspecies, nspec ", nspecies, NSPEC
     stop
  end if

  if (iH2 .eq. -100) then
     iH2 = get_species_index("H2")
     iO2 = get_species_index("O2")
     iN2 = get_species_index("N2")
     iH2O = get_species_index("H2O")
     iHO2 = get_species_index("HO2")
     iO = get_species_index("O")
     iH2O2 = get_species_index("H2O2")
     iH = get_species_index("H")
     iOH = get_species_index("OH")
  end if

  do i = state_l1, state_h1

     xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

     state(i,:) = 0.d0

  end do

end subroutine rns_initdata

