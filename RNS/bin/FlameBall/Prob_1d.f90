
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(1), probhi(1)

  integer untin,i

  namelist /fortin/ prob_type, pertmag, rfire, uinit, vinit, winit

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

!     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center(1) = 0.5d0*(problo(1) + probhi(1))

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

  use eos_module, only : eos_given_PTY
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC
  use chemistry_module, only : Patm, nspecies

  implicit none
  integer level, nscal
  integer lo(1), hi(1)
  integer state_l1,state_h1
  double precision state(state_l1:state_h1,NVAR)
  double precision time, delta(1)
  double precision xlo(1), xhi(1)
  
  integer :: i, n, iwrk
  double precision :: xcen, r, rfront
  double precision :: pmf_vals(NSPEC+3), Xt(nspec), Yt(nspec)
  double precision :: rhot, et, Pt, Tt, u1t, rwrk

  if (nspecies .ne. NSPEC) then
     write(6,*)"nspecies, nspec ", nspecies, NSPEC
     stop
  end if

  do i = state_l1, state_h1

     xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

     r = abs(xcen-center(1))
     rfront = rfire - r + 3.011d0 ! 3.011d0 is roughly the surface of fire for pmf.

     call pmf(rfront,rfront,pmf_vals,n)
     
     if (n .ne. nspec+3) then
        write(6,*)"n, nspec ", n, nspec
        stop
     end if

     Xt = pmf_vals(4:)

     Pt  = patm
     Tt  = pmf_vals(1)
     u1t = uinit
     
     call eos_given_PTY(rhot, et, Yt, Pt, Tt, Xt)

     state(i,URHO ) = rhot
     state(i,UMX  ) = rhot*u1t
     state(i,UEDEN) = rhot*(et + 0.5d0*u1t**2)
     state(i,UTEMP) = Tt

     do n=1,nspec
        state(i,UFS+n-1) = Yt(n)*rhot
     end do

  end do

end subroutine rns_initdata

