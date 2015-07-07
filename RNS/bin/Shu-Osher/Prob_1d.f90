
subroutine PROBINIT (init,name,namlen,problo,probhi)

  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(1), probhi(1)

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

  use eos_module, only : eos_get_T, gamma_const
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC

  implicit none
  integer level, nscal
  integer lo(1), hi(1)
  integer state_l1,state_h1
  double precision state(state_l1:state_h1,NVAR)
  double precision time, delta(1)
  double precision xlo(1), xhi(1)
  
  double precision xg, xcen, rho, u, p, Yt(2), Tt, et
  integer i, ii, n

  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  
  do i = state_l1, state_h1
     xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

     state(i,:) = 0.d0

     do ii = 1, 2
        xg = xcen + 0.5d0*delta(1)*gp(ii)

        if (xg < 1.d0) then
           rho = 3.857143d0
           u   = 2.629369d0
           p   = 10.33333d0
        else
           rho = 1.d0+0.2d0*sin(5.d0*xg)
           u = 0.d0
           p = 1.d0
        end if

        et = p/((gamma_const-1.d0)*rho)
        Yt = 0.5d0
        call eos_get_T(Tt, et, Yt)

        state(i,URHO ) = state(i,URHO ) + 0.5d0*rho
        state(i,UMX  ) = state(i,UMX  ) + 0.5d0*rho*u
        state(i,UEDEN) = state(i,UEDEN) + 0.5d0*rho*(et+0.5d0*u*u)
        state(i,UTEMP) = state(i,UTEMP) + 0.5d0*Tt
        do n=1, NSPEC
           state(i,UFS+n-1) = state(i,UFS+n-1) + 0.5d0*rho*Yt(n)
        end do
     end do            
  enddo

end subroutine rns_initdata


