
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ idir, max_denerr_lev, max_dengrad_lev, dengrad

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
  idir = 0
  max_denerr_lev = 10
  max_dengrad_lev = -1
  dengrad = 0.2

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center = 0.5d0*(problo + probhi)

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
  use meth_params_module, only : NVAR, URHO, UEDEN, UTEMP, UFS, NSPEC

  implicit none

  integer level, nscal
  integer lo(3), hi(3)
  integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision xlo(3), xhi(3), time, delta(3)
  double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
  
  ! local variables
  integer :: i, j, k, ii, jj, kk, n
  double precision :: xcen, ycen, zcen, xg, yg, zg, r
  double precision :: rhot, Pt, et, Tt, Yt(2)

  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: Pi = 3.1415926535897932d0

  if (NSPEC.ne. 2) then
     write(6,*)"nspec .ne. 2", NSPEC
     stop
  end if

  do k = lo(3), hi(3)
     zcen = xlo(3) + delta(3)*(dble(k-lo(3)) + 0.5d0)

     do j = lo(2), hi(2)
        ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)
        
        do i = lo(1), hi(1)
           xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)
           
           state(i,j,k,:) = 0.d0
           
           do kk = 1, 2
              zg = zcen + 0.5d0*delta(3)*gp(kk)
              do jj = 1, 2
                 yg = ycen + 0.5d0*delta(2)*gp(jj)
                 do ii = 1, 2
                    xg = xcen + 0.5d0*delta(1)*gp(ii)
                    
                    if (idir .eq. 1) then
                       r = sqrt((yg-center(2))**2 + (zg-center(3))**2)
                    else if (idir .eq. 2) then
                       r = sqrt((xg-center(1))**2 + (zg-center(3))**2)
                    else if (idir .eq. 3) then
                       r = sqrt((xg-center(1))**2 + (yg-center(2))**2)
                    else
                       r = sqrt((xg-center(1))**2 + (yg-center(2))**2 + (zg-center(3))**2)
                    end if
                    
                    if (r .gt. rpulse) then
                       rhot = rho0
                    else
                       rhot = rho0 + drho0*exp(-16.d0*r*r)*(cos(Pi*r))**6
                    end if
                    
                    Pt = (rhot/rho0)**gamma_const
                    
                    et = Pt/((gamma_const-1.d0)*rhot)
                    
                    Yt(1) = min(1.d0, rhot/(rho0+drho0))
                    Yt(2) = 1.d0 - Yt(1)
                    
                    call eos_get_T(Tt, et, Yt)
                    
                    state(i,j,k,URHO ) = state(i,j,k,URHO ) + 0.125d0*rhot
                    state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + 0.125d0*rhot*et
                    state(i,j,k,UTEMP) = state(i,j,k,UTEMP) + 0.125d0*Tt
                    do n=1, NSPEC
                       state(i,j,k,UFS+n-1) = state(i,j,k,UFS+n-1) + 0.125d0*rhot*Yt(n)
                    end do
                    
                 end do
              end do
           end do

        end do
     end do
  end do

end subroutine rns_initdata

