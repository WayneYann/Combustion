
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(1), probhi(1)

  center(1) = 0.5d0*(problo(1) + probhi(1))
  length(1) = probhi(1) - problo(1)

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
     state,state_l1,state_h1, &
     delta,xlo,xhi)

  use eos_module, only : gamma_const, eos_get_T
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC

  implicit none

  integer level, nscal
  integer lo(1), hi(1)
  integer state_l1,state_h1
  double precision xlo(1), xhi(1), time, delta(1)
  double precision state(state_l1:state_h1,NVAR)
  
  ! local variables
  integer :: i, ii, n
  double precision :: xcen, xg
  double precision :: rhot, Pt, et, Tt, Yt(2), ekt

  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: Pi = 3.1415926535897932d0

  if (NSPEC.ne. 2) then
     write(6,*)"nspec .ne. 2", NSPEC
     stop
  end if

  do i = lo(1), hi(1)
     xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

     state(i,:) = 0.d0

     do ii = 1, 2
        xg = xcen + 0.5d0*delta(1)*gp(ii)

        rhot = rho0 + drho0*sin(2.d0*Pi/length(1)*xg)

        Pt = p0
        
        et = Pt/((gamma_const-1.d0)*rhot)
        
        ekt = 0.5d0*u0*u0
        
        Yt(1) = 1.d0
        Yt(2) = 0.d0
        
        call eos_get_T(Tt, et, Yt)
        
        state(i,URHO ) = state(i,URHO ) + 0.5d0*rhot
        state(i,UMX  ) = state(i,UMX  ) + 0.5d0*rhot*u0
        state(i,UEDEN) = state(i,UEDEN) + 0.5d0*rhot*(et+ekt)
        state(i,UTEMP) = state(i,UTEMP) + 0.5d0*Tt
        do n=1, NSPEC
           state(i,UFS+n-1) = state(i,UFS+n-1) + 0.5d0*rhot*Yt(n)
        end do

     end do
  end do

end subroutine rns_initdata

