
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ prob_type, pertmag, rfire, uinit, vinit, winit, T0, T1, &
       max_denerr_lev, max_tracerr_lev, max_temperr_lev, temperr, tracerr, &
       fuel_name

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
  prob_type = 1

! problem type 1
  pertmag = 0.d0
  rfire   = 0.15d0
  uinit   = 0.d0
  vinit   = 0.d0
  winit   = 0.d0

  T0 = 1100.d0
  T1 = 1500.d0

  max_denerr_lev = -1
  max_temperr_lev = -1
  max_tracerr_lev = -1
  temperr = 1250.d0
  tracerr = 3.d-11

  fuel_name = "H2"

!     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center(1) = 0.5d0*(problo(1) + probhi(1))
  center(2) = 0.5d0*(problo(2) + probhi(2))

  Length(1) = probhi(1) - problo(1)
  Length(2) = probhi(2) - problo(2)

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

  use eos_module, only : eos_given_PTX
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UTEMP, UFS, NSPEC
  use chemistry_module, only : Patm, nspecies, get_species_index

  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)
  
  ! local variables
  integer :: i, j, n, ii, jj, nimages, iii, jjj, ifuel, iO2, iN2
  double precision :: xcen, ycen, xg, yg, r, rfront, xgi, ygi
  double precision :: pmf_vals(NSPEC+3), Xt(nspec), Yt(nspec)
  double precision :: rhot, et, Pt, Tt, u1t, u2t, kx, ky, Pi

!  integer, parameter :: ngp = 2
!  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
!  double precision, parameter :: wgp(2) = (/ 1.d0, 1.d0 /)
!
!  integer, parameter :: ngp = 3
!  double precision, parameter :: gp(3) = (/ -sqrt(0.6d0), 0.d0, sqrt(0.6d0) /)
!  double precision, parameter :: wgp(3) = (/ 5.d0/9.d0, 8.d0/9.d0, 5.d0/9.d0 /)  
!
  integer, parameter :: ngp = 4
  double precision, parameter :: gp(4) = (/  &
       -sqrt((3.d0+2.d0*sqrt(1.2d0))/7.d0), -sqrt((3.d0-2.d0*sqrt(1.2d0))/7.d0), &
       +sqrt((3.d0-2.d0*sqrt(1.2d0))/7.d0),  sqrt((3.d0+2.d0*sqrt(1.2d0))/7.d0) /)
  double precision, parameter :: wgp(4) = (/ &
       (18.d0-sqrt(30.d0))/36.d0, (18.d0+sqrt(30.d0))/36.d0, &
       (18.d0+sqrt(30.d0))/36.d0, (18.d0-sqrt(30.d0))/36.d0 /)

  double precision :: w

  if (nspecies .ne. NSPEC) then
     write(6,*)"nspecies, nspec ", nspecies, NSPEC
     stop
  end if

  ifuel = get_species_index(fuel_name)
  if (ifuel .le. 0) then
     stop "wrong fuel name"
  end if
  iO2 = get_species_index("O2")
  iN2 = get_species_index("N2")

  Pi = 4.d0*atan(1.d0)

  nimages = 0
  if (prob_type .eq. 4) nimages = 3

  do j = state_l2, state_h2
     ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)


     do i = state_l1, state_h1
        xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

        state(i,j,:) = 0.d0

        do jj = 1, ngp
           yg = ycen + 0.5d0*delta(2)*gp(jj)
           do ii = 1, ngp
              xg = xcen + 0.5d0*delta(1)*gp(ii)

              if (prob_type .eq. 1) then

                 r = sqrt((xg-center(1))**2 + (yg-center(2))**2)
                 rfront = rfire - r + 3.011d0 ! 3.011d0 is roughly the surface of fire for pmf.
              else if (prob_type .eq. 4) then

              else
                 write(6,*) "Unknown prob_type"
                 stop
              end if

              if (prob_type .eq. 1) then
                 call pmf(rfront,rfront,pmf_vals,n)
              
                 if (n .ne. nspec+3) then
                    write(6,*)"n, nspec ", n, nspec
                    stop
                 end if
              
                 Xt = pmf_vals(4:)
              end if
              
              if (prob_type .eq. 1) then
                 Pt  = patm
                 Tt  = pmf_vals(1)
                 u1t = uinit
                 u2t = vinit

              else if (prob_type .eq. 4) then
                 
                 Pt = Patm
                 Tt = T0

                 Xt = 0.0d0
                 Xt(ifuel) = 0.10d0
                 Xt(iO2) = 0.25d0

                 do jjj = -nimages, nimages
                    do iii = -nimages, nimages

                       xgi = xg + iii*Length(1)
                       ygi = yg + jjj*Length(2)
                       
                       r = sqrt(xgi**2+ygi**2)
                       
                       Pt = Pt    + 0.1d0*patm * exp(-(r / rfire)**2)
                       Tt = Tt       + (T1-T0) * exp(-(r / rfire)**2)
                       Xt(ifuel) = Xt(ifuel) + 0.025d0 * exp(-(r / rfire)**2)
                       Xt(iO2) = Xt(iO2) - 0.050d0 * exp(-(r / rfire)**2)
                       
                    end do
                 end do

                 kx = 2.d0*Pi/Length(1)
                 ky = 2.d0*Pi/Length(2)

                 u1t =  sin(kx*xg)*cos(ky*yg) * 300.d0
                 u2t = -cos(kx*xg)*sin(ky*yg) * 300.d0
              
                 Xt(iN2) = 1.0d0 - Xt(ifuel) - Xt(iO2)

              end if
              
              call eos_given_PTX(rhot, et, Yt, Pt, Tt, Xt)

              w = wgp(ii)*wgp(jj)*0.25d0

              state(i,j,URHO ) = state(i,j,URHO ) + w*rhot
              state(i,j,UMX  ) = state(i,j,UMX  ) + w*rhot*u1t
              state(i,j,UMY  ) = state(i,j,UMY  ) + w*rhot*u2t
              state(i,j,UEDEN) = state(i,j,UEDEN) + w*rhot*(et + 0.5d0*(u1t**2+u2t**2))
              state(i,j,UTEMP) = state(i,j,UTEMP) + w*Tt
              do n=1, NSPEC
                 state(i,j,UFS+n-1) = state(i,j,UFS+n-1) + w*rhot*Yt(n)
              end do
              
           end do
        end do

     end do
  end do
  
end subroutine rns_initdata

