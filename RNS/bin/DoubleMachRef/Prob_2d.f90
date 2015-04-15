
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  xshock0_hi = xshock0_lo + (probhi(2)-problo(2))/tan(thetashock)
  vshock_x = vshock/sin(thetashock)

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
  integer :: i, j, ii, jj
  double precision :: xcen, ycen, ei0, ei1, Y(2), xg, yg, w

  integer, parameter :: ngp = 2
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgp(2) = (/ 1.d0, 1.d0 /)

  logical, save :: first_call = .true.

  if (first_call) then

     first_call = .false.

     ei0 = p0/((gamma_const-1.d0)*rho0)
     ei1 = p1/((gamma_const-1.d0)*rho1)

     allocate(state0(NVAR))
     allocate(state1(NVAR))

     state0(URHO)  = rho0
     state0(UMX)   =  rho0*v0*sin(thetashock)
     state0(UMY)   = -rho0*v0*cos(thetashock)
     state0(UEDEN) = rho0*(ei0 + 0.5d0*v0**2)
     state0(UTEMP) = 0.d0
     Y(1) = 0.25d0
     Y(2) = 0.75d0
     call eos_get_T(state0(UTEMP), ei0, Y)
     state0(UFS)   = Y(1) * rho0
     state0(UFS+1) = Y(2) * rho0

     state1(URHO)  = rho1
     state1(UMX)   =  rho1*v1*sin(thetashock)
     state1(UMY)   = -rho1*v1*cos(thetashock)
     state1(UEDEN) = rho1*(ei1 + 0.5d0*v1**2)
     state1(UTEMP) = 0.d0
     Y(1) = 0.75d0
     Y(2) = 0.25d0
     call eos_get_T(state1(UTEMP), ei1, Y)
     state1(UFS)   = Y(1) * rho1
     state1(UFS+1) = Y(2) * rho1

  end if

  if (NSPEC.ne. 2) then
     write(6,*)"nspec .ne. 2", NSPEC
     stop
  end if

  do j = state_l2, state_h2
     ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)

     do i = state_l1, state_h1
        xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

        state(i,j,:) = 0.d0

        do jj = 1, ngp
           yg = ycen + 0.5d0*delta(2)*gp(jj)
           do ii = 1, ngp
              xg = xcen + 0.5d0*delta(2)*gp(ii)

              w = wgp(ii)*wgp(jj)*0.25d0

              if (yg .gt. tan(thetashock)*(xg-xshock0_lo)) then
                 state(i,j,:) = state(i,j,:) + w*state1
              else
                 state(i,j,:) = state(i,j,:) + w*state0
              end if
           end do
        end do

     end do
  end do

end subroutine rns_initdata

