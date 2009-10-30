
c-----------------------------------------------------------------------
c     -----------------------------------------------------------
c     This routine is called at problem initialization time
c     and when restarting from a checkpoint file.
c     The purpose is (1) to specify the initial time value
c     (not all problems start at time=0.0) and (2) to read
c     problem specific data from a namelist or other input
c     files and possibly store them or derived information
c     in FORTRAN common blocks for later use.
c     
c     
c     INPUTS/OUTPUTS:
c     
c     init      => TRUE if called at start of problem run
c     FALSE if called from restart
c     strttime <=  start problem with this time variable
c     
c     -----------------------------------------------------------
      subroutine PROBINIT (problo,probhi)
      implicit none

      integer untin
      DOUBLE PRECISION problo, probhi
      DOUBLE PRECISION rn,permin,permax,xtmp,pert

      include 'nums.fi'
      include 'probdata.fi'
      include 'sndata.fi'
      include 'nkbrn.fi'
      include 'bc.fi'

      namelist /fortin/ probtype,nspec,nreac,forceInFlow,splitx,
     &                  xfrontw,dpdt_factor,
     &                  max_temp_lev,tempgrad,
     &                  max_spec_lev,specgrad,
     &                  flameloc,pertamp,nfreq

      integer lo, hi, i
      integer maxlen
      integer ierr
      parameter (maxlen=32)

      integer MAXPHASE
      parameter (MAXPHASE=10000)

      DOUBLE PRECISION xfrontw_frac
      integer n

c     Default values

      probtype = 1

      forceInFlow = .FALSE.
      vorterr = 1.0d20
      flametracval = 0.0001d0
      specgrad = 50.0D0
      max_spec_lev = 0
      tempgrad = 50.0D0
      max_temp_lev = 0
      max_vort_lev = 0
      splitx = 0.0D0
      flameloc = 0.d0
      nfreq = 0
      pertamp = 0.d0
      dpdt_factor = 0.0D0

      xfrontw_frac = 1.0D0/8.0D0
      xfrontw = (probhi-problo)*xfrontw_frac

      untin = 9
      open(untin,file='probin',form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      xfrontw = xfrontw / 4.d0

      Density = 1
      Temp = 2
      RhoH = 3
      FirstSpec = 4
      LastSpec = FirstSpec + (nspec-1)

c     Load domain dimensions into common, and set up boundary functions
      domnlo = problo
      domnhi = probhi

c     Pi = 4.0D0*atan2(1.d0,1.d0)
c     call blutilinitrand(111397)
c     do i = 1, nfreq
c        call blutilrand(rn)
c        ranampl(i,1)=2.d0*(rn-0.5d0)
c        call blutilrand(rn)
c        ranphse(i,1,1)=2.d0*Pi*rn
c     end do
c     permin = 100.d0
c     permax = -100.d0
c     do i = 0, MAXPHASE
c        xtmp = 2.d0*Pi*dble(i)/dble(MAXPHASE)
c        pert = 0.d0
c        do n=1,nfreq
c           pert = pert+sin(2.0D0*Pi*
c    &           dble(n)*xtmp+ranphse(n,1,1))
c    &           *ranampl(n,1)
c        end do
c        permin = min(permin,pert)
c        permax = max(permax,pert)
c     end do
c     pertamp = 2.d0*pertamp/(permax-permin)
C     
      call setupbc()
      bcinit = .true.


      end
c-----------------------------------------------------------------------
      subroutine setupbc()
      implicit none
      include 'nums.fi'
      include 'nkbrn.fi'
      include 'bc.fi'
      include 'probdata.fi'
      include 'sndata.fi'
      integer lo, hi
      data lo / 1 /
      data hi / 1 /
      integer n, zone, ierr, cnt, nspec1
      integer NDUM
      parameter (NDUM = 27 + MAXSPEC )
      DOUBLE PRECISION state1d(NDUM)
      DOUBLE PRECISION Y(Nzones*MAXSPEC)
      DOUBLE PRECISION T(Nzones)
      DOUBLE PRECISION R(Nzones)
      DOUBLE PRECISION U(Nzones)
      DOUBLE PRECISION V(Nzones)
      integer i
c
      if ( probtype .eq. 1 ) then

         Y(1) = 0.5
         Y(2) = 0.
         Y(3) = 0.5
         Y(4) = 0.
         Y(5) = 1.0
         Y(6) = 0.

         R(1) = 5.0e7
         R(2) = -1.

         T(1) = 1.0e7
         T(2) = 4.0e9

         U(1) = 9.11e4
         U(2) = 0.
C     
         do zone = 1, Nzones
            do n = 1, Nspec
               Y_bc(n, zone) = Y(n + (zone-1)*Nspec)
            end do
            T_bc(zone) = T(zone)
            u_bc(zone) = U(zone)
            v_bc(zone) = V(zone)
         end do
C     
C     Zone 2:
C     

         if ( R(2).gt. 0.d0 ) then
            call nb_pgivenrty(
     &           Pamb,1,1,
     &           R(2),1,1,
     &           T(2),1,1,
     &           Y(1+(2-1)*Nspec),1,1,
     &           lo, hi)
            R_bc(2) = R(2)
!     
!     Zone 1:
!     
            R_bc(1) = R(2)      ! initial guess
            call nb_rgivenpty(
     &           R_bc(1),1,1,
     &           Pamb,1,1,
     &           T(1),1,1,
     &           Y(1+(1-1)*Nspec),1,1,
     &           lo, hi, cnt)

         else


            call nb_pgivenrty(
     &           Pamb,1,1,
     &           R(1),1,1,
     &           T(1),1,1,
     &           Y(1+(1-1)*Nspec),1,1,
     &           lo, hi)
            R_bc(1) = R(1)
!     
!     Zone 1:
!     
            R_bc(2) = R(1)      ! initial guess
            call nb_rgivenpty(
     &           R_bc(2),1,1,
     &           Pamb,1,1,
     &           T(2),1,1,
     &           Y(1+(2-1)*Nspec),1,1,
     &           lo, hi, cnt)
         end if

!         print *, 'R_bc = ', R_bc
!         print *, 'Pamb = ', Pamb
!         print *, 'cnt = ', cnt

         do zone = 1, Nzones
            call nb_hgivenrty(
     &           h_bc(zone),1,1,
     &           R_bc(zone),   1,1,
     &           T_bc(zone),   1,1,
     &           Y(1+(zone-1)*Nspec),1,1,
     &           lo,hi)
         end do

      end if

      end
c-----------------------------------------------------------------------
      integer function getZone(x)
      implicit none
      DOUBLE PRECISION x
      include 'nkbrn.fi'
      include 'bc.fi'
      include 'probdata.fi'

      if ( probtype .eq. 1) then
         getZone = 1
         if ( x .GE. splitx ) getZone = 2
      else if ( probtype .eq. 2 ) then
         if ( x.LE.0.d0 ) then
            getZone=1
         else
            getZone=2
         end if
      else
         call bl_error('INVALID PROBTYPE')
      end if
      end
c-----------------------------------------------------------------------
      subroutine bcfunction(x, time, u, rho, Yl, T, h)
      implicit none
      DOUBLE PRECISION x, time, u, rho, Yl(*), T, h
      include 'nums.fi'
      include 'nkbrn.fi'
      include 'sndata.fi'
      include 'bc.fi'
      include 'probdata.fi'
      DOUBLE PRECISION s, eta, mdotL, mdotR, mdotb
      DOUBLE PRECISION zblend
      external zblend
      integer n, getZone, zone, lo, hi, zL, zR
c     integer Nspec, cnt
      integer cnt
      data lo / 1 /
      data hi / 1 /
c     A handy statement function: here blend goes from 0 to 1 in x at
c     splitx, over a width of xfrontw.  
c     Some boundary conditions may not be not time-dependent, but rather
c     just constant values for each orientation.  We can precompute
c     then (in setupbc), here just retrieve the appropriate values.
c     However, we do not know the size of the domain here, so we cannot
c     call setupbc values ourselves. Need to assume someone else did it.
      if (.not. bcinit) then
         call bl_error('Need to initialize boundary condition function')
      end if

      if ( probtype.eq.1 ) then
         zL = getZone(domnlo)
         zR = getZone(domnhi)
         eta = zblend(x)
         T = (1.0D0-eta)*T_bc(zL) + eta*T_bc(zR)
         do n = 1,Nspec
            Yl(n) = (1.0D0-eta)*Y_bc(n,zL) +eta*Y_bc(n,zR)
         end do
c     initial guess
         rho = (1.0D0-eta)*R_bc(zL)+eta*R_bc(zR)
         call nb_rgivenpty(
     &        rho,1,1,
     &        Pamb,1,1,
     &        T,1,1,
     &        Yl,1,1,
     &        lo, hi, cnt)
         call nb_hgivenrty(
     &        h,1,1,
     &        rho,1,1,
     &        T,1,1,
     &        Yl,1,1,
     &        lo,hi)
         mdotL = u_bc(zL)*rho
         mdotR = u_bc(zR)*rho
         mdotB = (1.0D0-eta)*mdotL + eta*mdotR
         u = mdotb/rho
      else if ( probtype .eq. 2 ) then
         zone = getZone(x)
         if ( zone .eq. 1 ) then
            u=u_bc(zone)
            rho=R_bc(zone)
            T=T_bc(zone)
            h=h_bc(zone)
            do n = 1, Nspec
               Yl(n)= y_bc(n,zone)
            end do
         else
            call bl_error(
     &           "called bcfunction and got zone=2 and protype=2"
     &           )
         end if
      else
         call bl_error('INVALID PROBTYPE')
      end if
      end
c----------------------------------------------------------------------
c     A handy statement function: here blend goes from 0 to 1 in x at
c     splitx, over a width of xfrontw.  
      DOUBLE PRECISION function zblend(x)
      DOUBLE PRECISION x
      include 'probdata.fi'
      zblend = 0.5D0*(1.0D0+TANH((x-splitx)/(xfrontw)))
      end
c---------------------------------------------------------------------
c     -----------------------------------------------------------
c     This routine is called at problem setup time and is used
c     to initialize data on each grid.  The velocity field you
c     provide does not have to be divergence free and the pressure
c     field need not be set.  A subsequent projection iteration
c     will define aa divergence free velocity field along with a
c     consistant pressure.
c     
c     NOTE:  all arrays have 1.0D0 cell of ghost zones surrounding
c     the grid interior.  Values in these cells need not
c     be set here.
c     
c     INPUTS/OUTPUTS:
c     
c     time      => time at which to init data             
c     lo,hi     => index limits of grid interior (cell centered)
c     nscal     => number of scalar quantities.  You should know
c     this already!
c     vel      <=  Velocity array
c     scal     <=  Scalar array
c     press    <=  Pressure array
c     dx        => cell size
c     right hand corner of grid.  (does not include
c     ghost region).
c     -----------------------------------------------------------
      subroutine initdata(nx,vel,scal,Ydot,dx)

      implicit none

      integer              nx
      DOUBLE PRECISION     dx
      DOUBLE PRECISION    vel(-1:nx)
      DOUBLE PRECISION   scal(-1:nx  ,nscal)
      DOUBLE PRECISION   Ydot(0 :nx-1,nspec)

      include 'nums.fi'
      include 'nkbrn.fi'
      include 'sndata.fi'
      include 'bc.fi'
      include 'probdata.fi'

      DOUBLE PRECISION y, x, eta,yloloc,yhiloc,pert
      integer zL, zR 
      integer NDUM
      parameter (NDUM = 27 + MAXSPEC)
      DOUBLE PRECISION mdotl, mdotr, mdotb
      DOUBLE PRECISION state1d(NDUM)
      integer k, j, i, n, cnt
      DOUBLE PRECISION  u, v, rho, Yl(maxspec), T, h
      DOUBLE PRECISION  time

      time = 0.0

      do n = 1,nspec
      do i = 0,nx-1
        Ydot(i,n) = 0.d0
      end do
      end do

      do i = 0,nx-1
            x = (dble(i)+0.5D0)*dx
            call bcfunction(x, time, u, rho, Yl, T, h)
            scal(i,Temp) = T
            do n = 1, Nspec
               scal(i,FirstSpec+n-1) = rho*Yl(n)
            end do
            scal(i,Density) = rho
            scal(i,RhoH) = rho*h
            vel(i) = u
      end do

      do n = 1,nscal
        scal(-1,n) = scal(0   ,n)
        scal(nx,n) = scal(nx-1,n)
      end do

      vel(-1) = vel(0)
      vel(nx) = vel(nx-1)

      end

c     subroutine bl_error(str)
c     character(*) str
c     print *,str
c     stop
c     end
c     
c     subroutine bl_pd_is_ioproc(bl)
c     integer bl
c     bl = 1
c     end

      subroutine nb_get_maxvode(imvs)
      imvs = 1000
      end

