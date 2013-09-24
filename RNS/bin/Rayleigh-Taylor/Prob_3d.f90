
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ prob_type, frac, rho_1, rho_2, p0_base, pertmag, &
       dengrad, max_dengrad_lev

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)
  
  if (namlen .gt. maxlen) then
     call bl_error('probin file name too long')
  end if
  
  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
  
  ! set namelist defaults here
  prob_type = 0
  frac = 0.5d0
  rho_1 = 1.0d0
  rho_2 = 2.0d0
  p0_base = 5.0d0
  pertmag = 0.0d0
  
  dengrad = 0.01
  max_dengrad_lev = 5

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! set local variable defaults
  center(1) = frac*(problo(1)+probhi(1))
  center(2) = frac*(problo(2)+probhi(2))
  center(3) = frac*(problo(3)+probhi(3))
  
  L_x = probhi(1) - problo(1)
  L_y = probhi(2) - problo(2)

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

  use eos_module, only : gamma_const, eos_get_T
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, NSPEC

  implicit none

  integer level, nscal
  integer lo(3), hi(3)
  integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision xlo(3), xhi(3), time, delta(3)
  double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
  
  ! local variables
  integer :: i, j, k
  double precision :: xcen,ycen,r2d,zcen,ei,pres,presmid,pertheight, Y(2)
  double precision :: vx, vy, vz, ek
  double precision, parameter :: ZERO=0.d0, HALF=0.5d0, &
       PI = 3.141592653589793238462643383279502884197d0

  double precision :: r 
  integer :: nseed
  integer, allocatable :: seed(:)
  logical, save :: first_call = .true.

  if (first_call) then
     call random_seed(SIZE=nseed)
     allocate(seed(nseed))
     call random_seed(GET=seed)
     seed = seed + sum(lo) + sum(hi)
     call random_seed(PUT=seed)
     deallocate(seed)
     first_call = .false.
  end if

  if (NSPEC.ne. 2) then
     write(6,*)"nspec .ne. 2", NSPEC
     stop
  end if

  if (prob_type .eq. 0) then

     presmid  = p0_base - rho_1*center(3)

     state(:,:,:,UMX)   = ZERO
     state(:,:,:,UMY)   = ZERO
     state(:,:,:,UMZ)   = ZERO

     do k = state_l3, state_h3
        zcen = (k+HALF)*delta(3)

        do j = state_l2, state_h2           
           do i = state_l1, state_h1
              if (zcen .lt. center(3)) then
                 pres = p0_base - rho_1*zcen
                 state(i,j,k,UEDEN) = pres / (gamma_const - 1.0d0)
              else
                 pres = presmid - rho_2*(zcen-center(3))
                 state(i,j,k,UEDEN) = pres / (gamma_const - 1.0d0)
              end if
           end do
        end do
     end do

     do k = state_l3, state_h3
        zcen = (k+HALF)*delta(3)

        do j = state_l2, state_h2
           ycen = (j+HALF)*delta(2)

           if (pertmag .gt. 0.0d0) then
              call random_number(r)
              ycen = ycen + pertmag*(r-0.5d0)*delta(2)
           end if

           do i = state_l1, state_h1
              xcen = (i+HALF)*delta(1)

              if (pertmag .gt. 0.0d0) then
                 call random_number(r)
                 xcen = xcen + pertmag*(r-0.5d0)*delta(1)
              end if

              r2d = min(sqrt((xcen-center(1))**2+(ycen-center(2))**2), 0.5d0*L_x)
              pertheight = 0.5d0 - 0.01d0*cos(2.0d0*PI*r2d/L_x)
              state(i,j,k,URHO) = rho_1 + ((rho_2-rho_1)/2.0d0)* &
                   (1.d0+tanh((zcen-pertheight)/0.005d0))

              Y(1) = (rho_2 - state(i,j,k,URHO)) / (rho_2-rho_1)
              Y(2) = 1.d0 - Y(1)
              ei = state(i,j,k,UEDEN) / state(i,j,k,URHO)
              state(i,j,k,UTEMP) = ZERO
              call eos_get_T(state(i,j,k,UTEMP), ei, Y)

              state(i,j,k,UFS  ) = state(i,j,k,URHO) * Y(1)
              state(i,j,k,UFS+1) = state(i,j,k,URHO) * Y(2)

           enddo
        enddo
     enddo

  else if (prob_type .eq. 1) then

     do k = state_l3, state_h3
        zcen = xlo(3) + (dble(k-lo(3))+HALF)*delta(3)

        do j = state_l2, state_h2
           ycen = xlo(2) + (dble(j-lo(2))+HALF)*delta(2)
           
           do i = state_l1, state_h1
              xcen = xlo(1) + (dble(i-lo(1))+HALF)*delta(1)

              if (zcen .lt. 0.d0) then
                 state(i,j,k,URHO) = rho_1
                 pres = P0_base - 0.1d0*rho_1*ycen
                 Y(1) = 1.0d0
                 Y(2) = 0.0d0
              else
                 state(i,j,k,URHO) = rho_2
                 pres = P0_base - 0.1d0*rho_2*ycen
                 Y(1) = 1.0d0
                 Y(2) = 0.0d0
              end if

              ei = pres / ((gamma_const - 1.d0) * state(i,j,k,URHO))
              vx = 0.d0
              vy = 0.d0
              vz = 0.01d0*(1.d0+cos(4.d0*PI*xcen))*(1.d0+cos(4.d0*PI*ycen))* &
                   (1.d0+cos(3.d0*PI*zcen))/4.d0
              ek = 0.5d0*(vx**2 + vy**2 + vz**2)

              state(i,j,k,UMX  ) = state(i,j,k,URHO) * vx
              state(i,j,k,UMY  ) = state(i,j,k,URHO) * vy
              state(i,j,k,UMZ  ) = state(i,j,k,URHO) * vz
              state(i,j,k,UEDEN) = state(i,j,k,URHO) * (ei+ek)

              state(i,j,k,UTEMP) = ZERO
              call eos_get_T(state(i,j,k,UTEMP), ei, Y)

              state(i,j,k,UFS  ) = state(i,j,k,URHO) * Y(1)
              state(i,j,k,UFS+1) = state(i,j,k,URHO) * Y(2)

           end do
        end do
     end do

  else

     write(6,*)"Unknown prob_type", prob_type
     stop

  end if

end subroutine rns_initdata

