
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ prob_type, turbfile, pamb, phi_in, T_in, vn_in, T_co, vn_co, &
       splitx, xfrontw, Tfrontw, blobr, blobx, bloby, blobT, inflow_period, inflow_vnmag, &
       splity, yfrontw, turb_boost_factor, &
       max_tracerr_lev, tracerr, max_vorterr_lev, vorterr, max_tempgrad_lev, tempgrad

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

  prob_type = 2

  pamb = 4.053d7  ! 40 Patm
  
  phi_in = 0.2d0   ! mole fraction of CH3OCH3 in fuel
  T_in   = 400.d0  ! temperature of fuel
  vn_in  = 5.12d3  ! fuel injection velocity
  
  T_co  = 1525.d0   ! temperature of air
  vn_co = 5.12d2    ! air injection velocity

  splitx  = 0.00569d0  ! where fuel and air split
  xfrontw = .2d0       ! controls the width of split
  Tfrontw = 0.0075d0

  inflow_period = 1.d-5 ! period of sinusoidal variation of inflow velocity
  inflow_vnmag  = 1.d3  ! magnitude of the variation
  
  blobr = -1.d0
  blobx = 0.d0
  bloby = 0.027d0
  blobT = 1500.d0

  splity  = 0.001d0
  yfrontw = 0.0004d0

  turb_boost_factor = 1.d0

  max_tracerr_lev = -1
  tracerr = 1.d-8

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
  integer :: i, j, k, n, ii, jj, kk
  double precision :: xcen, ycen, zcen, xg, yg, zg, r
  double precision Yt(NSPEC)
  double precision rhot,u1t,u2t,u3t,Tt,et, eta
  integer :: iwrk
  double precision :: rwrk
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgt = 0.125d0

  if (nspecies .ne. NSPEC) then
     write(6,*)"nspecies, nspec ", nspecies, NSPEC
     stop
  end if

  if (.not. dmejet_initialized) then
     call init_DME_jet()
  end if

  do k = state_l3, state_h3
     zcen = xlo(3) + delta(3)*(dble(k-lo(3)) + 0.5d0)

     do j = state_l2, state_h2
        ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)
        
        do i = state_l1, state_h1
           xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)
           
           state(i,j,k,:) = 0.d0

           do kk = 1, 2
              zg = zcen + 0.5d0*delta(3)*gp(kk)
              do jj = 1, 2
                 yg = ycen + 0.5d0*delta(2)*gp(jj)
                 do ii = 1, 2
                    xg = xcen + 0.5d0*delta(1)*gp(ii)
                    
                    r = sqrt(xg*xg+yg*yg)
                    
                    ! prob_type 2 only
                    
                    eta = 0.5d0*(1.d0-tanh((r -splitx)/xfrontw))
                    if ((zg-splity) < 5.d0*yfrontw) then
                       eta = eta * 0.5d0*(1.d0-tanh((zg -splity)/yfrontw))
                    else
                       eta = 0.d0
                    end if
                    
                    do n=1,nspecies
                       Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
                    end do
                    Tt  = eta * T_in + (1.d0-eta) * T_co
                    u1t = 0.d0
                    u2t = 0.d0
                    u3t = eta * vn_in + (1.d0-eta) * vn_co
                    

                    CALL CKRHOY(pamb,Tt,Yt,IWRK,RWRK,rhot)
                    call CKUBMS(Tt,Yt,IWRK,RWRK,et)
                    
                    state(i,j,k,URHO ) = state(i,j,k,URHO ) + wgt*rhot
                    state(i,j,k,UMX  ) = state(i,j,k,UMX  ) + wgt*rhot*u1t
                    state(i,j,k,UMY  ) = state(i,j,k,UMY  ) + wgt*rhot*u2t
                    state(i,j,k,UMZ  ) = state(i,j,k,UMZ  ) + wgt*rhot*u3t
                    state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + wgt*rhot*(et + 0.5d0*(u1t**2+u2t**2+u3t**2))
                    state(i,j,k,UTEMP) = state(i,j,k,UTEMP) + wgt*Tt
                    do n=1, NSPEC
                       state(i,j,k,UFS+n-1) = state(i,j,k,UFS+n-1) + wgt*rhot*Yt(n)
                    end do
                    
                 end do
              end do
           end do
           
        end do
     end do
  end do

end subroutine rns_initdata

