
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use eos_module, only : gamma_const
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(1), probhi(1)

  integer untin,i

  namelist /fortin/ prob_type, xsep, p_l, u_l, rho_l, p_r, u_r, rho_r

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

  prob_type = 0
  xsep = 0.5d0

  p_l = 1.0               ! left pressure (erg/cc)
  u_l = 0.0               ! left velocity (cm/s)
  rho_l = 1.0             ! left density (g/cc)

  p_r = 0.1               ! right pressure (erg/cc)
  u_r = 0.0               ! right velocity (cm/s)
  rho_r = 0.125           ! right density (g/cc)

!     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center(1) = 0.5d0*(problo(1)+probhi(1))

  select case (prob_type) 
     case (1)
        rho_l = 1.d0
        u_l   = 0.75d0
        p_l   = 1.d0
        rho_r = 0.125
        u_r   = 0.d0
        p_r   = 0.1d0
        xsep  = 0.3
     case (2)
        rho_l = 1.d0
        u_l   = -2.d0
        p_l   = 0.4d0
        rho_r = 1.d0
        u_r   = 2.d0
        p_r   = 0.4d0
        xsep  = 0.5d0
     case (3)
        rho_l = 1.d0
        u_l   = 1.d0
        p_l   = 1.d-6
        rho_r = 1.d0
        u_r   = -1.d0
        p_r   = 1.d-6
        xsep  = 0.5d0
     case (30)    ! case 3a
        rho_l = 1.d0
        u_l   = -19.59745d0
        p_l   = 1000.d0
        rho_r = 1.d0
        u_r   = -19.59745d0
        p_r   = 0.01d0
        xsep  = 0.8d0
     case (4)
        rho_l = 5.99924d0
        u_l   = 19.5975d0
        p_l   = 460.894d0
        rho_r = 5.99242d0
        u_r   = -6.19633d0
        p_r   = 46.095d0
        xsep  = 0.4d0
     case (5)
        rho_l = 1.4d0
        u_l   = 0.d0
        p_l   = 1.d0
        rho_r = 1.d0
        u_r   = 0.d0
        p_r   = 1.d0
        xsep  = 0.5d0
     case (6)
        rho_l = 1.4d0
        u_l   = 0.1d0
        p_l   = 1.d0
        rho_r = 1.d0
        u_r   = 0.1d0
        p_r   = 1.d0
        xsep  = 0.5d0
     case (7)   ! case peak
        rho_l = 0.1261192d0
        u_l   = 8.9047029d0
        p_l   = 782.92899d0
        rho_r = 6.591493d0
        u_r   = 2.2654207d0
        p_r   = 3.154487d0
        xsep  = 0.5d0
  end select

!     compute the internal energy (erg/cc) for the left and right state
  rhoe_l = p_l/(gamma_const - 1.d0)
  rhoe_r = p_r/(gamma_const - 1.d0)

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

  use eos_module, only : eos_get_T
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC

  implicit none
  integer level, nscal
  integer lo(1), hi(1)
  integer state_l1,state_h1
  double precision state(state_l1:state_h1,NVAR)
  double precision time, delta(1)
  double precision xlo(1), xhi(1)
  
  double precision xcen, T, e, Y(NSPEC)
  integer i
  
  do i = state_l1, state_h1
     xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)
            
     if (xcen <= xsep) then
        state(i,URHO ) = rho_l
        state(i,UMX  ) = rho_l*u_l
        state(i,UEDEN) = rhoe_l + 0.5d0*rho_l*u_l*u_l
        e = rhoe_l/rho_l

        Y(1) = 6.d0/7.d0
        Y(2) = 1.d0/7.d0
        call eos_get_T(T, e, Y)
        state(i,UTEMP) = T
        state(i,UFS:UFS+1) = state(i,URHO) * Y

     else
        state(i,URHO ) = rho_r
        state(i,UMX  ) = rho_r*u_r
        state(i,UEDEN) = rhoe_r + 0.5d0*rho_r*u_r*u_r
        e = rhoe_r/rho_r

        Y(1) = 1.d0/7.d0
        Y(2) = 6.d0/7.d0
        call eos_get_T(T, e, Y)
        state(i,UTEMP) = T
        state(i,UFS:UFS+1) = state(i,URHO) * Y
     endif
  enddo

end subroutine rns_initdata


