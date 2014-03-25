
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(1), probhi(1)

  integer untin,i

  namelist /fortin/ prob_type, pertmag, rfire, uinit, vinit, winit, T0, T1, &
       max_denerr_lev, max_tracerr_lev, tracerr

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
  prob_type = 1

! problem type 1
  pertmag = 0.d0
  rfire   = 0.15d0
  uinit   = 0.d0
  vinit   = 0.d0
  winit   = 0.d0

  T0 = 1100.d0
  T1 = 1500.d0

  max_denerr_lev = -1
  max_tracerr_lev = -1
  tracerr = 3.d-11

!     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center(1) = 0.5d0*(problo(1) + probhi(1))

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
  
  integer :: i, n, ii, nimages, iii, iH2, iO2, iN2
  double precision :: xcen, xg, r, rfront, xgi
  double precision :: pmf_vals(NSPEC+3), Xt(nspec), Yt(nspec)
  double precision :: rhot, et, Pt, Tt, u1t

!  integer, parameter :: ngp = 2
!  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
!  double precision, parameter :: wgp(2) = (/ 1.d0, 1.d0 /)
!
!  integer, parameter :: ngp = 3
!  double precision, parameter :: gp(3) = (/ -sqrt(0.6d0), 0.d0, sqrt(0.6d0) /)
!  double precision, parameter :: wgp(3) = (/ 5.d0/9.d0, 8.d0/9.d0, 5.d0/9.d0 /)  
!
  integer, parameter :: ngp = 4
  double precision, parameter :: gp(4) = (/  &
       -sqrt((3.d0+2.d0*sqrt(1.2d0))/7.d0), -sqrt((3.d0-2.d0*sqrt(1.2d0))/7.d0), &
       +sqrt((3.d0-2.d0*sqrt(1.2d0))/7.d0),  sqrt((3.d0+2.d0*sqrt(1.2d0))/7.d0) /)
  double precision, parameter :: wgp(4) = (/ &
       (18.d0-sqrt(30.d0))/36.d0, (18.d0+sqrt(30.d0))/36.d0, &
       (18.d0+sqrt(30.d0))/36.d0, (18.d0-sqrt(30.d0))/36.d0 /)

  double precision :: w

  if (nspecies .ne. NSPEC) then
     write(6,*)"nspecies, nspec ", nspecies, NSPEC
     stop
  end if

  iH2 = get_species_index("H2")
  iO2 = get_species_index("O2")
  iN2 = get_species_index("N2")

  nimages = 0
  if (prob_type .eq. 4) nimages = 3

  do i = state_l1, state_h1

     xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

     state(i,:) = 0.d0

     do ii = 1, ngp
        xg = xcen + 0.5d0*delta(1)*gp(ii)

        if (prob_type  .eq. 1) then

           r = abs(xg-center(1))
           rfront = rfire - r + 3.011d0 ! 3.011d0 is roughly the surface of fire for pmf.

        else if (prob_type .eq. 4) then

        else
           write(6,*) "Unknown prob_type"
           stop
        end if
        
        if (prob_type .eq. 1) then
           call pmf(rfront,rfront,pmf_vals,n)
              
           if (n .ne. nspec+3) then
              write(6,*)"n, nspec ", n, nspec
              stop
           end if
           
           Xt = pmf_vals(4:)
        end if

        if (prob_type .eq. 1) then

           Pt  = patm
           Tt  = pmf_vals(1)
           u1t = uinit
           
        else if (prob_type .eq. 4) then
                 
           Pt = Patm
           Tt = T0

           Xt = 0.0d0
           Xt(iH2) = 0.10d0
           Xt(iO2) = 0.25d0
           
           do iii = -nimages, nimages

              xgi = xg + iii*Length(1)
                       
              r = abs(xgi)
                       
              Pt = Pt    + 0.1d0*patm * exp(-(r / rfire)**2)
              Tt = Tt       + (T1-T0) * exp(-(r / rfire)**2)
              Xt(iH2) = Xt(iH2) + 0.025d0 * exp(-(r / rfire)**2)
              Xt(iO2) = Xt(iO2) - 0.050d0 * exp(-(r / rfire)**2)
                       
           end do

           u1t = 0.0
              
           Xt(iN2) = 1.0d0 - Xt(iH2) - Xt(iO2)

        end if
     
        call eos_given_PTX(rhot, et, Yt, Pt, Tt, Xt)

        w = wgp(ii)*0.5d0

        state(i,URHO ) = state(i,URHO ) + w*rhot
        state(i,UMX  ) = state(i,UMX  ) + w*rhot*u1t
        state(i,UEDEN) = state(i,UEDEN) + w*rhot*(et + 0.5d0*(u1t**2))
        state(i,UTEMP) = state(i,UTEMP) + w*Tt
        do n=1, NSPEC
           state(i,UFS+n-1) = state(i,UFS+n-1) + w*rhot*Yt(n)
        end do

     end do
  end do

end subroutine rns_initdata

