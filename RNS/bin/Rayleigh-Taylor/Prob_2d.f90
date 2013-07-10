
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

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
  
  L_x = probhi(1) - problo(1)

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
     state,state_l1,state_l2,state_h1,state_h2, &
     delta,xlo,xhi)

  use eos_module, only : gamma_const, eos_get_T
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UTEMP, UFS, NSPEC

  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)
  
  ! local variables
  integer :: i, j, ncell_x
  double precision :: xcen,ycen,ei,pres,presmid,pertheight, Y(2)
  double precision :: vx, vy, ek
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
     seed = seed + lo(1) + lo(2) + hi(1) + hi(2)
     call random_seed(PUT=seed)
     deallocate(seed)
     first_call = .false.
  end if

  if (NSPEC.ne. 2) then
     write(6,*)"nspec .ne. 2", NSPEC
     stop
  end if

  if (prob_type .eq. 0) then

     ncell_x = nint(L_x/delta(1))

     presmid  = p0_base - rho_1*center(2)

     state(:,:,UMX)   = ZERO
     state(:,:,UMY)   = ZERO

     do j = state_l2, state_h2
        ycen = (j+HALF)*delta(2)

        do i = state_l1, state_h1

           if (ycen .lt. center(2)) then
              pres = p0_base - rho_1*ycen
              state(i,j,UEDEN) = pres / (gamma_const - 1.0d0)
           else
              pres = presmid - rho_2*(ycen-center(2))
              state(i,j,UEDEN) = pres / (gamma_const - 1.0d0)
           end if

        end do
     end do

     do j = state_l2, state_h2
        ycen = (j+HALF)*delta(2)

        do i = state_l1, state_h1

           if (i .lt. ncell_x/2) then
              xcen = (i+HALF)*delta(1)
           else
              xcen = ((ncell_x-i-1)+HALF)*delta(1)
           end if

           ! we explicitly make the perturbation symmetric here
           ! -- this prevents the RT from bending.
           pertheight = 0.01d0*HALF*(cos(2.0d0*PI*xcen/L_x) + &
                cos(2.0d0*PI*(L_x-xcen)/L_x)) + 0.5d0

           if (pertmag .gt. 0.0d0) then
              call random_number(r)
              pertheight = pertheight * (1.d0 + r*1.d-10)
           end if

           state(i,j,URHO) = rho_1 + ((rho_2-rho_1)/2.0d0)* &
                (1+tanh((ycen-pertheight)/0.005d0))

           Y(1) = (rho_2 - state(i,j,URHO)) / (rho_2-rho_1)
           Y(2) = 1.d0 - Y(1)
           ei = state(i,j,UEDEN) / state(i,j,URHO)
           state(i,j,UTEMP) = ZERO
           call eos_get_T(state(i,j,UTEMP), ei, Y)

           state(i,j,UFS  ) = state(i,j,URHO) * Y(1)
           state(i,j,UFS+1) = state(i,j,URHO) * Y(2)

        enddo
     enddo

  else if (prob_type .eq. 1) then

     do j = state_l2, state_h2
        ycen = xlo(2) + (dble(j-lo(2))+HALF)*delta(2)

        do i = state_l1, state_h1
           xcen = xlo(1) + (dble(i-lo(1))+HALF)*delta(1)

           if (ycen .lt. 0.d0) then
              state(i,j,URHO) = rho_1
              pres = P0_base - 0.1d0*rho_1*ycen
              Y(1) = 1.0d0
              Y(2) = 0.0d0
           else
              state(i,j,URHO) = rho_2
              pres = P0_base - 0.1d0*rho_2*ycen
              Y(1) = 1.0d0
              Y(2) = 0.0d0
           end if

           ei = pres / ((gamma_const - 1.d0) * state(i,j,URHO))
           vx = 0.d0
           vy = 0.01d0*(1.d0+cos(4.d0*PI*xcen))*(1.d0+cos(3.d0*PI*ycen))/4.d0
           ek = 0.5d0*(vx**2 + vy**2)

           state(i,j,UMX  ) = state(i,j,URHO) * vx
           state(i,j,UMY  ) = state(i,j,URHO) * vy
           state(i,j,UEDEN) = state(i,j,URHO) * (ei+ek)

           state(i,j,UTEMP) = ZERO
           call eos_get_T(state(i,j,UTEMP), ei, Y)

           state(i,j,UFS  ) = state(i,j,URHO) * Y(1)
           state(i,j,UFS+1) = state(i,j,URHO) * Y(2)

        end do
     end do

  else

     write(6,*)"Unknown prob_type", prob_type
     stop

  end if

end subroutine rns_initdata


! ::: -----------------------------------------------------------

subroutine rns_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc)
  
  use meth_params_module, only : NVAR
  use probdata_module
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)
  
  integer i, j, n
  
  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,n), &
          adv_l1,adv_l2,adv_h1,adv_h2, &
             domlo,domhi,delta,xlo,bc(1,1,n))
  enddo
  
  do n=1,NVAR
     !        XLO
     if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
        call bl_error('SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) ')
     end if
     
     !        XHI
     if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
        call bl_error('SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) ')
     end if
         
     !        YLO
     if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
        call bl_error('SHOULD NEVER GET HERE bc(2,1,n) .eq. EXT_DIR) ')
     end if
         
     !        YHI
     if ( bc(2,2,n).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
        call bl_error('SHOULD NEVER GET HERE bc(2,2,n) .eq. EXT_DIR) ')
     end if

  end do

end subroutine rns_hypfill

! ::: -----------------------------------------------------------

subroutine rns_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)
      
  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     print *,'SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     print *,'SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_denfill
