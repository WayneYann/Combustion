
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


  iH2 = get_species_index("H2")
  iO2 = get_species_index("O2")
  iN2 = get_species_index("N2")
  iH2O = get_species_index("H2O")
  iHO2 = get_species_index("HO2")
  iO = get_species_index("O")
  iH2O2 = get_species_index("H2O2")
  iH = get_species_index("H")
  iOH = get_species_index("OH")
  
  allocate(X_reac(nspecies))
  allocate(X_prod(nspecies))
  allocate(X_intm(nspecies))

  X_reac = 0.d0
  X_reac(iH2) = 0.296d0
  X_reac(iO2) = 0.148d0
  X_reac(iN2) = 0.556d0

  X_prod = 0.d0
  X_prod(iH2O) = 0.296d0
  X_prod(iO2)  = 0.000d0
  X_prod(iN2)  = 0.556d0

  X_intm = 0.d0
  X_intm(iHO2)  = 0.0001d0
  X_intm(iO)    = 0.0001d0
  X_intm(iH2O2) = 0.0001d0
  X_intm(iH)    = 0.0100d0
  X_intm(iOH)   = 0.0100d0

  T_reac = 298.d0
  T_prod = 2300.d0
  massFlux = 0.07d0
  pres = 1.d0 * patm

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

  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC
  use chemistry_module, only : nspecies, get_species_index
  use inflow_module

  implicit none
  integer level, nscal
  integer lo(1), hi(1)
  integer state_l1,state_h1
  double precision state(state_l1:state_h1,NVAR)
  double precision time, delta(1)
  double precision xlo(1), xhi(1)
  
  integer :: i, n, iwrk, iN2
  double precision :: Xt(nspec), Yt(nspec), Yti(nspecies), rwrk
  double precision :: xcen, rhot,Tt,et,Pt, eta, eta01, eta00

  if (nspecies .ne. NSPEC) then
     write(6,*)"nspecies, nspec ", nspecies, NSPEC
     stop
  end if

  if (.not. allocated(inflow_state)) call init_inflow()

  iN2 = get_species_index("N2")

  do i = state_l1, state_h1

     xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0) - init_flm_position

     state(i,:) = 0.d0

     Pt = pres

     if (xcen .lt. -10.d0*init_flm_width) then
        Tt = T_reac
        Xt = X_reac
        Xt = Xt/sum(Xt)
        CALL CKXTY (Xt, IWRK, RWRK, Yt)
     else if (xcen .gt. 10.d0*init_flm_width) then
        Tt = T_prod
        Xt = X_prod
        Xt = Xt/sum(Xt)
        CALL CKXTY (Xt, IWRK, RWRK, Yt)
     else
        eta = tanh(xcen/init_flm_width)
        eta01 = (eta+1.d0)*0.5d0
        Tt = (1.d0-eta01)*T_reac + eta01*T_prod
        Xt = (1.d0-eta01)*X_reac + eta01*X_prod
        Xt = Xt/sum(Xt)
        CALL CKXTY (Xt, IWRK, RWRK, Yti)
        ! add flame
        eta00 = 1.d0 - abs(eta)
        Xt = Xt + eta00*X_intm
        Xt = Xt/sum(Xt)
        CALL CKXTY (Xt, IWRK, RWRK, Yt)
        Yt = Yt*(1.d0-Yti(iN2))/(1.d0-Yt(iN2))
        Yt(iN2) = Yti(iN2)          
     end if
     
     CALL CKRHOY(Pt,Tt,Yt,IWRK,RWRK,rhot)
     call CKUBMS(Tt,Yt,IWRK,RWRK,et)
     
     state(i,URHO) = rhot
     state(i,UMX)  = massFlux
     state(i,UEDEN) = rhot*et + 0.5d0*state(i,UMX)**2/state(i,URHO)
     state(i,UTEMP) = Tt
     do n=1,nspecies
        state(i,UFS-1+n) = Yt(n)*rhot
     end do

  end do

end subroutine rns_initdata

