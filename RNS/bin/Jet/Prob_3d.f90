
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ turbfile, pamb, phi_in, T_in, vn_in, T_co, vn_co, &
       splitx, xfrontw, turb_boost_factor, &
       max_vorterr_lev, vorterr, max_tempgrad_lev, tempgrad

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

  pamb = 4.053d7  ! 40 Patm
  
  phi_in = 0.2d0   ! mole fraction of CH3OCH3 in fuel
  T_in   = 400.d0  ! temperature of fuel
  vn_in  = 5.12d3  ! fuel injection velocity
  
  T_co  = 1525.d0   ! temperature of air
  vn_co = 5.12d2    ! air injection velocity

  splitx  = 0.00569d0  ! where fuel and air split
  xfrontw = .2d0       ! controls the width of split

  turb_boost_factor = 1.d0

  max_vorterr_lev = -1
  vorterr = 5.d4

  max_tempgrad_lev = -1
  tempgrad = 10.d0

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center = 0.5d0*(problo+probhi)

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
     state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
     delta,xlo,xhi)

  use eos_module, only : eos_given_PTX
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, NSPEC
  use chemistry_module, only : Patm, nspecies

  implicit none

  integer level, nscal
  integer lo(3), hi(3)
  integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision xlo(3), xhi(3), time, delta(3)
  double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
  
  ! local variables
  integer :: i, j, k, n
  double precision Yt(NSPEC)
  double precision rhot,u1t,u2t,u3t,Tt,et
  integer :: iwrk
  double precision :: rwrk

  if (nspecies .ne. NSPEC) then
     write(6,*)"nspecies, nspec ", nspecies, NSPEC
     stop
  end if

  if (.not. jet_initialized) then
     call init_jet()
  end if

  !$omp parallel do private(i,j,k,n,Yt,rhot,u1t,u2t,u3t,Tt,et,iwrk,rwrk) collapse(2)
  do k = state_l3, state_h3
     do j = state_l2, state_h2        
        do i = state_l1, state_h1

           do n=1,nspecies
              Yt(n) = air_Y(n)
           end do
           Tt = T_co
           u1t = 0.d0
           u2t = 0.d0
           u3t = vn_co

           CALL CKRHOY(pamb,Tt,Yt,IWRK,RWRK,rhot)
           call CKUBMS(Tt,Yt,IWRK,RWRK,et)
                    
           state(i,j,k,URHO ) = rhot
           state(i,j,k,UMX  ) = rhot*u1t
           state(i,j,k,UMY  ) = rhot*u2t
           state(i,j,k,UMZ  ) = rhot*u3t
           state(i,j,k,UEDEN) = rhot*(et + 0.5d0*(u1t**2+u2t**2+u3t**2))
           state(i,j,k,UTEMP) = Tt
           do n=1, NSPEC
              state(i,j,k,UFS+n-1) = rhot*Yt(n)
           end do
           
        end do
     end do
  end do
  !$omp end parallel do

end subroutine rns_initdata

