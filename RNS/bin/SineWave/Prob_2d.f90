
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  center(1:2) = 0.5d0*(problo + probhi)
  length(1:2) = probhi - problo

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
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC

  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)
  
  ! local variables
  integer :: i, j, ii, jj, n
  double precision :: xcen, ycen, xg, yg
  double precision :: rhot, Pt, et, Tt, Yt(2), ekt

  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: Pi = 3.1415926535897932d0

  if (NSPEC.ne. 2) then
     write(6,*)"nspec .ne. 2", NSPEC
     stop
  end if

  do j = lo(2), hi(2)
     ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)

     do i = lo(1), hi(1)
        xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)
        
        state(i,j,:) = 0.d0

        do jj = 1, 2
           yg = ycen + 0.5d0*delta(2)*gp(jj)
           do ii = 1, 2
              xg = xcen + 0.5d0*delta(1)*gp(ii)

              rhot = rho0 + drho0*sin(2.d0*Pi/length(1)*xg)

              Pt = p0

              et = Pt/((gamma_const-1.d0)*rhot)

              ekt = 0.5d0*u0*u0

              Yt(1) = 1.d0
              Yt(2) = 0.d0

              call eos_get_T(Tt, et, Yt)
              
              state(i,j,URHO ) = state(i,j,URHO ) + 0.25d0*rhot
              state(i,j,UMX  ) = state(i,j,UMX  ) + 0.25d0*rhot*u0
              state(i,j,UEDEN) = state(i,j,UEDEN) + 0.25d0*rhot*(et+ekt)
              state(i,j,UTEMP) = state(i,j,UTEMP) + 0.25d0*Tt
              do n=1, NSPEC
                 state(i,j,UFS+n-1) = state(i,j,UFS+n-1) + 0.25d0*rhot*Yt(n)
              end do

           end do
        end do

     end do
  end do

end subroutine rns_initdata

