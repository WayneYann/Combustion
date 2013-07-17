
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ xsep, rho0, rho1, pl, pr

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
  xsep = 0.1d0
  rho0 = 1.0d0
  rho1 = 1.0d0
  pl   = 1000.d0
  pr   = 0.01d0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center(1) = 0.5d0*(problo(1)+probhi(1))
  center(2) = 0.5d0*(problo(2)+probhi(2))

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
  integer :: i, j, ipert, jpert
  double precision :: xcen, ycen, pt, et, Y(2), rhot

  if (NSPEC.ne. 2) then
     write(6,*)"nspec .ne. 2", NSPEC
     stop
  end if

  ipert = int(2.d0*xsep/delta(1))
  jpert = int(center(2)/delta(2))

  do j = state_l2, state_h2
     ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)

     do i = state_l1, state_h1
        xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

        if (xcen < xsep) then
           pt = pl
        else
           pt = pr
        end if

        if (i .eq. ipert .and. j .eq. jpert) then
           rhot = rho1
        else
           rhot = rho0
        end if

        et = pt / ((gamma_const-1.d0)*rhot)

        state(i,j,URHO ) = rhot
        state(i,j,UMX  ) = 0.d0
        state(i,j,UMY  ) = 0.d0
        state(i,j,UEDEN) = rhot*et
        state(i,j,UTEMP) = 0.d0
        Y(1) = 0.9d0
        Y(2) = 1.d0 - Y(1)
        call eos_get_T(state(i,j,UTEMP), et, Y)
        state(i,j,UFS)   = Y(1) * rhot
        state(i,j,UFS+1) = Y(2) * rhot

     end do
  end do

end subroutine rns_initdata

