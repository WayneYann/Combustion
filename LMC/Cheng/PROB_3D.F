#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROB_AMR_F.H"
#include "PROB_F.H"
#include "ArrayLim.H"
#include "ChemKinDriver_F.H"
#include "HEATTRANSFER_F.H"
#include "infl_frc.H"

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

#define SDIM 3
#define USE_PMF_AVG
c#define SWIRL

      subroutine FORT_HACK(lo, hi, dat, DIMS(dat),
     &                     datc, DIMS(datc), nc, mapO2N, lenmap, ratio)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(dat)
      integer DIMDEC(datc)
      integer nc, lenmap, mapO2N(lenmap), ratio
      REAL_T dat(DIMV(dat),nc)
      REAL_T datc(DIMV(datc),nc)
      call bl_abort('Should not be in FORT_HACK for this problem...')
      end

      subroutine FORT_SET_PROB_SPEC(fuel, oxid, prod, numspec)
      implicit none
#include "probdata.H"
      integer fuel, oxid, prod, numspec
      fuelID = fuel + 1
      oxidID = oxid + 1
      prodID = prod + 1
      end

c ::: -----------------------------------------------------------
c ::: This routine is called at problem initialization time
c ::: and when restarting from a checkpoint file.
c ::: The purpose is (1) to specify the initial time value
c ::: (not all problems start at time=0.0) and (2) to read
c ::: problem specific data from a namelist or other input
c ::: files and possibly store them or derived information
c ::: in FORTRAN common blocks for later use.
c ::: 
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: init      => TRUE if called at start of problem run
c :::              FALSE if called from restart
c ::: strttime <=  start problem with this time variable
c ::: 
c ::: -----------------------------------------------------------
      subroutine FORT_PROBINIT (init,name,namlen,problo,probhi)
      implicit none
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(SDIM), probhi(SDIM)

#include "probdata.H"
#include "htdata.H"
#include "ckdwrk.H"
#include "bc.H"
#include "INFL_FORCE_F.H"

      integer lo(SDIM), hi(SDIM), i, j
      INTEGER dimFile(3)
      REAL_T dxFile(3)
      REAL_T rand

      namelist /fortin/ vorterr, temperr, adverr, tempgrad, 
     &                  flametracval, probtype,
     &		          max_temp_lev, max_vort_lev, 
     &                  traceSpecVal,
     &                  max_nozzle_lev, refine_nozzle, refine_nozzle_x, 
     &                  refine_nozzle_y, refine_nozzle_z, nozzle_width,
     &                  max_stick_lev, refine_stick, refine_stick_x, 
     &                  refine_stick_z, direction,
     &                  forceInflow, numInflPlanesStore, strmwse_dir,
     $                  forceLo, forceHi, flct_file, turb_scale,Vin, rhot, zhot,
     $                  zstandoff, Thot, swK, swW, stTh, Ro, Rf, wallTh,
     &                  zBL, dBL, stBL, hotspot_r, hotspot_z, hotspot_dz,
     &                  hotspot_T, hotspot_htc, hotspot_w, vheight,
     &                  Vco_l,Vco_r,tVco_l,tVco_r,anisotsc
      namelist /heattransin/ pamb, dpdt_factor
c
c      Build `probin' filename -- the name of file containing fortin namelist.
c
      integer maxlen
      parameter (maxlen=256)

      character probin*(maxlen)
      REAL_T y, yEval, pmf_vals(maxspec+3)
      integer n

      if (namlen .gt. maxlen) then
         call bl_abort('probin file name too long')
      end if

      if (namlen .eq. 0) then
         namlen = 6
         probin(1:namlen) = 'probin'
      else
         do i = 1, namlen
            probin(i:i) = char(name(i))
         end do
      endif

      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      
      anisotsc = 1.3d0
      vorterr = 1.e20
      temperr = zero
      adverr = 1.e20
      tempgrad  = 50.0d0
      flametracval = 0.0001d0
      probtype = 1
      max_temp_lev = 0
      max_vort_lev = 0
      traceSpecVal = 1.d-10
      refine_nozzle = 0
      refine_nozzle_z = zero
      refine_nozzle_y = zero
      refine_nozzle_x = zero
      nozzle_width = -one
      max_nozzle_lev = 100
      refine_stick = 0
      refine_stick_x = zero
      refine_stick_z = zero
      max_stick_lev = 100
      pamb = 101325.d0
      dpdt_factor = 0.3d0
      direction = "z"
      vheight = zero

      add_turb = .FALSE.
      forceInflow = .FALSE.
      numInflPlanesStore = -1
      forceLo = .TRUE.
      forceHi = .FALSE.
      strmwse_dir = FLCT_ZVEL
      flct_file = ""
      turb_scale = 1

c     Note: for setup with no coflow, set Ro=Rf+wallth
      Rf = .025d0
      wallTh = .001d0
      stTh = .002d0
      Ro = .05d0
      dBL = .01d0
      zBL = tenth*Rf
      stBL = half*dBL
      Vin = 5.d0

      Vco_l  = 0.1d0
      Vco_r  = 0.1d0
      tVco_l = 0.0d0
      tVco_r = 0.0d0

      Thot = 2.d3
      rhot = Rf
      zhot = Rf
      zstandoff = zero
      swK = zero
      swW = 0.6d0

      read(untin,fortin)
c      write(6,fortin)
      
      if (max_vort_lev.lt.0) max_vort_lev=max_temp_lev
      
      read(untin,heattransin)
c      write(6,heattransin)

   10 close(unit=untin)

      if (forceInflow .eqv. .FALSE.) then
         forceLo = .FALSE.
         forceHi = .FALSE.
      else
#if defined(BL_DO_FLCT)
         if (flct_file.ne."") then
#define FF_UNIT 20
            write(6,*) '...initializing turbulence, reading header info'
            open(FF_UNIT, file=flct_file, form='unformatted')
            call RD_FLCTHD(FF_UNIT, dimFile, probSizeFile, dxFile)
            close(FF_UNIT)
         endif
#endif
      endif
      convVel = Vin

c     Load domain dimensions into common, and set up boundary functions
      domnlo(1) = problo(1)
      domnlo(2) = problo(2)
      domnlo(3) = problo(3)
      domnhi(1) = probhi(1)
      domnhi(2) = probhi(2)
      domnhi(3) = probhi(3)
      
      call setupbc()
      bcinit = .true.
      
      end

c ::: -----------------------------------------------------------

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
#define ZLO 4
#define ZHI 5

      subroutine setupbc()
      implicit none
#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "htdata.H"
      
      REAL_T Patm, x, pmf_vals(maxspec+3), sum, Y_pmf(maxspec)
      integer i, zone, n, lo(SDIM), hi(SDIM), iN2, getZone
      character*(maxspnml) name
      data lo / 1, 1, 1 /
      data hi / 1, 1, 1 /
      
c     Zone definitions

#define BL_FUELPIPE 1
#define BL_COFLOW   2
#define BL_STICK    3
#define BL_WALL     4
#define BL_AMBIENT  5
#define BL_VOLUME   6
#define BL_PIPEEND  7

      REAL_T X_CH4(Nzones),X_O2(Nzones),X_N2(Nzones),T(Nzones)
      REAL_T Xt(maxspec)
c     ZONES:            1                   2         3             4           5     6         7
#ifdef SWIRL
c     phi = 0.8
      data X_CH4 /  0.0774907749077d0,   zero, 0.0774907749077d0,   zero,    zero,   zero ,   zero /
      data X_O2  /   0.193726937269d0,  0.21d0, 0.193726937269d0,  0.21d0,  0.21d0, 0.21d0,  0.21d0 /
      data T     /        3.d2,          3.d2,         3.d2,        3.d2,    3.d2,   3.d2 ,   3.d2 /
#else
c     phi = 1.0
c     data X_CH4 /  0.095023d0,   zero,  0.095023d0,   zero,   zero,   zero  , zero   /
c     data X_O2  /  0.190045d0, 0.21d0,  0.190045d0,  0.21d0,  0.21d0, 0.21d0, 0.21d0 /
c     phi = 0.6
c     data X_CH4 /  0.0592662276576d0,   zero, 0.0592662276576d0,   zero,    zero,   zero , zero  /
c     data X_O2  /   0.197554092192d0,  0.21d0, 0.197554092192d0,  0.21d0,  0.21d0, 0.21d0, 0.21d0/
c     phi = 0.7
      data X_CH4 /  0.0684676292501d0,   zero, 0.0684676292501d0,   zero,    zero,   zero , zero   /
      data X_O2  /   0.195621797857d0,  0.21d0, 0.195621797857d0,  0.21d0,  0.21d0, 0.21d0, 0.21d0 /
c     phi = 0.75
c     data X_CH4 /  0.0730011587486d0,   zero, 0.0730011587486d0,   zero,    zero,   zero , zero   /
c     data X_O2  /   0.194669756663d0,  0.21d0, 0.194669756663d0,  0.21d0,  0.21d0, 0.21d0, 0.21d0 /
      data T     /        3.d2,          3.d2,         3.d2,        3.d2,    3.d2,   3.d2,   3.d2  /
#endif
      Patm = pamb / 101325.0d0
      iN2 = -1
      do n = 1,Nspec
         call get_spec_name(name, n)
         if (name .eq. 'N2' ) iN2 = n
      end do
         
c     Initialize zone values
      do zone = 1,Nzones
         sum = zero
         do n = 1,Nspec
            Xt(n) = zero
            call get_spec_name(name, n)
            if (name .EQ. 'CH4') Xt(n) = X_CH4(zone)
            if (name .EQ. 'O2')  Xt(n) = X_O2(zone)
            if (name .NE. 'N2')  sum = sum + Xt(n)
         end do
         Xt(iN2) = one - sum
         CALL CKXTY (Xt, IWRK, RWRK, Y_bc(0,zone))
         T_bc(zone) = T(zone)
         call FORT_RHOfromPTY(lo,hi,
     &        rho_bc(zone),ARLIM(lo),ARLIM(hi),
     &        T_bc(zone),ARLIM(lo),ARLIM(hi),
     &        Y_bc(0,zone),ARLIM(lo),ARLIM(hi), Patm)
         call FORT_HMIXfromTY(lo,hi,
     &        h_bc(zone),ARLIM(lo),ARLIM(hi),
     &        T_bc(zone),ARLIM(lo),ARLIM(hi),
     &        Y_bc(0,zone),ARLIM(lo),ARLIM(hi))
      end do

      end
         
c ::: -----------------------------------------------------------

      integer function getZone(x, y, z)
      implicit none
#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
      REAL_T x, y, z, r, x1
c     ZONE DEF:
c     BL_FUELPIPE, BL_COFLOW, BL_STICK, BL_WALL, BL_AMBIENT, BL_VOLUME
      getZone = BL_VOLUME
      if ( ((direction.eq."z") .and. (z.le.domnlo(3))).or.
     &     ((direction.eq."y") .and. (y.le.domnlo(2))).or.
     &     ((direction.eq."x") .and. (x.le.domnlo(1)))) then

         if (direction.eq."z") then
            r = SQRT(x*x + y*y)
            x1 = x
         else if (direction.eq."y") then
            r = SQRT(x*x + z*z)
            x1 = x
         else
            r = SQRT(y*y + z*z)
            x1 = y
         end if

         if (r.le.Rf) then
            getZone = BL_FUELPIPE         
            if (ABS(x1).lt.half*stTh) getZone = BL_STICK
         else if (r.le.Rf+wallTh) then
            getZone = BL_PIPEEND
         else if (r.le.Ro) then
            getZone = BL_COFLOW
         else if (r.ge.Ro+wallTh) then
            getZone = BL_AMBIENT
         else
            getZone = BL_WALL
         end if

      end if
      end
      
c
c     Compute the average velocity in a dx(1)*dx(2) square
c     x,y  position on entry plane
c     Rf     fuel tube radius
c     wallTh fuel tube thickness
c     Ro     oxidizer tube radius
c     stTh   thickness of stick crossing origin
c     ypls   dimensionless laminar sublayer thickness
c     ypt    dimensionless height of transition layer
c     nu     viscosity of fluid
c     F      flux of flux (cm**2/s)
c     r0     radius where log-law is replaced by parabola
c
      subroutine VAvg(V,x,y,dx,Rf,wallTh,Ro,stTh,dBL,stBL,
     &     Vf,Vco_l,Vco_r,tVco_l,tVco_r,time,swK,swW)
      implicit none
      REAL_T V(SDIM)
      REAL_T x,y,dx(SDIM),Rf,wallTh,Ro,stTh,dBL,stBL
      REAL_T Vf,Vco_l,Vco_r,tVco_l,tVco_r,time,swK,swW

      REAL_T incx,incy,x1,x2,eta
      REAL_T Vincx,Vincy,Vincz,Vswirl,Vin,r
      integer i,j,M
      parameter (M=10)

c     size of subvolumes, there are M*M of them
      incx = dx(1) / M
      incy = dx(2) / M
      V(1) = 0.d0
      V(2) = 0.d0
      V(3) = 0.d0
      do i=1,M
         do j=1,M
            x1 = x - 0.5d0*dx(1) + (i-0.5d0)*incx
            x2 = y - 0.5d0*dx(2) + (j-0.5d0)*incy
            r = sqrt(x1*x1 + x2*x2)
            Vincx = 0.d0
            Vincy = 0.d0
            Vincz = 0.d0
            if (r.le.Ro) then
#ifdef SWIRL
c               eta = r/Ro
c               Vincz = Vf*(0.63336d0 - 2.3529d0*eta + 25.595d0*eta**2 - 115.84d0*eta**3
c     $              + 255.47d0*eta**4 - 260.06d0*eta**5+ 99.786d0*eta**6)

#else
               Vincz = Vin*TANH(two*(one - r/Ro)/dBL)
#endif
            else if (r.ge.Ro+wallTh) then
               if (time .le. tVco_l) then
                  Vincz = Vco_l
               else if (time .ge. tVco_r) then
                  Vincz = Vco_r
               else
                  Vincz = Vco_l+(time-tVco_l)*(Vco_r-Vco_l)/(tVco_r-tVco_l)
               endif
            endif

            if (stTh.gt.zero) then
#ifdef SWIRL
               call bl_abort('better not be here')
#endif
               if (ABS(x1).le.half*stTh.and.r.le.Ro+wallTh) then
                  Vincz = 0.d0
               elseif(r.le.Ro+wallTh)then
                  Vincz = Vincz * TANH(2.d0*(ABS(x1)-half*stTh)/stBL)
               end if
            end if

#ifdef SWIRL
            if (r.le.Ro) then
c
c     JBB fit
c              Vswirl = swK*1.1d0*EXP(-20*(.95-r/Rf)**2)*TANH(40*(1.d0-r/Rf))
c
c     Original swirl parameterization
c              eta = r/Ro
c              if ((r.gt.Rf*(one-swW)).and.(r.le.Rf)) then
c                 Vswirl = swK*(SIN(Pi*(one-r/Rf)/swW))**2
c                 Vincx = Vswirl*(-x2/r)
c                 Vincy = Vswirl*( x1/r)
c              end if
c
c     Fit from averaged EBHC calc
c              Vswirl = Vin*(-0.067641d0 + 2.0113d0*eta - 19.305d0*eta**2 + 72.404d0*eta**3
c     &                     - 126.14d0*eta**4 + 111.66d0*eta**5 - 31.098d0*eta**6)
c              Vradial =Vin*(-0.0022569d0 + 0.3562d0*eta - 0.84751d0*eta**2 + 1.8533d0*eta**3
c     &                    - 5.1925d0*eta**4 + 7.7047d0*eta**5 - 3.8714d0*eta**6)
c              Vincx = Vswirl*(-x2/r) + Vradial*x1/r
c              Vincy = Vswirl*( x1/r) + Vradial*x2/r
            end if
#endif
              
            V(1) = V(1) + Vincx
            V(2) = V(2) + Vincy
            V(3) = V(3) + Vincz
         end do
      end do
      V(1) = V(1) / (M*M)
      V(2) = V(2) / (M*M)
      V(3) = V(3) / (M*M)
      end

c ::: -----------------------------------------------------------
      
      subroutine bcfunction(orient,x,y,z,time,u,v,w,rho,Yl,T,h,dx,getuvw)

      implicit none

      integer orient
      REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h, dx(SDIM)
      logical getuvw

#include "htdata.H"
#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"

      REAL_T Vbndry(3), eta, Patm, r, rn
      integer n, getZone, zone, lo(SDIM), hi(SDIM)
      data lo / 1, 1, 1 /
      data hi / 1, 1, 1 /

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if
      
      if (orient .le. ZHI) then

         zone = getZone(x,y,z)
         rho = rho_bc(zone)
         do n = 0, Nspec-1
            Yl(n) = Y_bc(n,zone)
         end do
         T = T_bc(zone)
         h = h_bc(zone)
         Patm = pamb / 101325.0d0

         if (zone.eq.BL_FUELPIPE) then

            r = SQRT(x**2 + y**2)
#ifdef SWIRL
            rn = r/Ro

c     Fit as of 4/29/03
            eta = 9.920628558770E-001
     &           -1.065813294748E-001 *rn**1
     &           +1.800416951718E+000 *rn**2
     &           -1.158582069687E+001 *rn**3
     &           +3.588317172715E+001 *rn**4
     &           -5.816747423739E+001 *rn**5
     &           +4.632391073835E+001 *rn**6
     &           -1.450956235808E+001 *rn**7
            
c     Ramp down the outter value to zero, 6/2/03
            eta = eta*0.5*(1-TANH((rn-0.92)/.03))

            
c            eta = 1.002d0 - 0.098868d0*r + 0.84254d0*r**2 - 2.6231d0*r**3 +
c     $           3.0126d0*r**4 - 1.9876d0*r**5 + 0.52246d0*r**6
c
            eta = min(max(eta,0.d0),1.d0)
#else
            if (stTh.gt.0.d0) then
               eta = 1.d0
            else
               eta = 0.5d0*(1 - TANH((r - 0.02d0)/.002d0))
            endif
#endif
            do n = 0, Nspec-1
               Yl(n) = (1.d0-eta)*Y_bc(n,BL_COFLOW) + eta*Y_bc(n,BL_FUELPIPE)
            end do

            call FORT_RHOfromPTY(lo,hi,
     &           rho,ARLIM(lo),ARLIM(hi),
     &           T,ARLIM(lo),ARLIM(hi),
     &           Yl,ARLIM(lo),ARLIM(hi), Patm)
            call FORT_HMIXfromTY(lo,hi,
     &           h,ARLIM(lo),ARLIM(hi),
     &           T,ARLIM(lo),ARLIM(hi),
     &           Yl,ARLIM(lo),ARLIM(hi))
         endif

         if (getuvw) then
            u = zero
            v = zero
            w = zero
            call VAvg(Vbndry,x,y,dx,Rf,wallTh,Ro,stTh,dBL,stBL,
     &           Vin,Vco_l,Vco_r,tVco_l,tVco_r,time,swK,swW)
            u = Vbndry(1)
            v = Vbndry(2)
            w = Vbndry(3)
         endif
      else
         write(6,*) 'No boundary condition for orientation = ', orient
         call bl_abort(" ")
      endif
      end

      REAL_T function shapet(r)
      REAL_T r
      if(r.le.1.d0)then
         shapet = tanh(100.d0*(1-r))
      else
         shapet = 0.d0
      endif
      end

      REAL_T function shapen(r)
      REAL_T r
      if(r.le.1.d0)then
         shapen = tanh(100.d0*(1-r))
      else
         shapen = 0.d0
      endif
      end

      REAL_T function turbSclT(r)
      REAL_T r, eta
#include"probdata.H"
      turbSclT = 0.d0
      eta = r/Ro
      if(r.lt.Ro)then
#ifdef SWIRL
         turbSclT = Vin*(0.10278d0 + 0.34958d0 *eta - 0.86045d0 *eta**2
     %     - 6.9187d0 *eta**3 + 30.903d0 *eta**4 - 31.34d0 *eta**5 + 8.588d0 *eta**6)

#else
         turbSclT = 1.d0
#endif
      endif
      end

      REAL_T function turbSclR(r)
      REAL_T r, eta
#include"probdata.H"
      turbSclR = 0.d0
      eta = r/Ro
      if(r.lt.Ro)then
#ifdef SWIRL
         turbSclR=Vin*(0.086871d0 + 1.5445d0*eta - 17.316d0*eta**2 + 79.469d0*eta**3
     &      - 164.43d0*eta**4 + 160.85d0*eta**5 - 60.233d0*eta**6)

#else
         turbSclR = (0.25d0 + 0.95d0*tanh((1.d0-r/Ro)/.0502362048))/1.2d0
#endif
      endif
      end

      REAL_T function turbSclX(x,y,z)
#include "probdata.H"
      REAL_T x, y, z, shapet, eta, Rmid
      eta = SQRT(x*x+y*y)/Ro
      turbSclX = shapet(eta)
      if (stBL.gt.zero) then
         turbSclX = turbSclX*TANH(4.d0*MAX(0.d0,ABS(x)-0.5d0*stTh)/stBL)
      end if
      end

      REAL_T function turbSclY(x,y,z)
      REAL_T x, y, z, turbSclX
      turbSclY = turbSclX(x,y,z)
      end

      REAL_T function turbSclZ(x,y,z)
#include "probdata.H"
      REAL_T x, y, z, shapen, eta, Rmid
      turbSclZ = zero
      eta = SQRT(x*x+y*y)/Ro
#ifdef SWIRL
      turbSclZ = one
      return
c      if(eta.lt.one)then
c         turbSclZ = Vin*(-0.020318d0 + 2.18d0*eta - 18.157d0*eta**2 + 66.538d0*eta**3 - 
c     &         112.55d0*eta**4 + 94.931d0*eta**5 - 32.329d0*eta**6)
c      endif
#else
      if(eta.lt.one)then
         turbSclZ = 2.d0 -1.d0 * tanh((1-eta)/.0428222042d0)
      endif
      if (stBL.gt.zero) then
         turbSclZ = turbSclZ*TANH(4.d0*MAX(0.d0,ABS(x)-0.5d0*stTh)/stBL)
      end if
#endif
      end

c ::: -----------------------------------------------------------
c ::: This routine is called at problem setup time and is used
c ::: to initialize data on each grid.  The velocity field you
c ::: provide does not have to be divergence free and the pressure
c ::: field need not be set.  A subsequent projection iteration
c ::: will define aa divergence free velocity field along with a
c ::: consistant pressure.
c ::: 
c ::: NOTE:  all arrays have one cell of ghost zones surrounding
c :::        the grid interior.  Values in these cells need not
c :::        be set here.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: level     => amr level of grid
c ::: time      => time at which to init data             
c ::: lo,hi     => index limits of grid interior (cell centered)
c ::: nscal     => number of scalar quantities.  You should know
c :::		   this already!
c ::: vel      <=  Velocity array
c ::: scal     <=  Scalar array
c ::: press    <=  Pressure array
c ::: delta     => cell size
c ::: xlo,xhi   => physical locations of lower left and upper
c :::              right hand corner of grid.  (does not include
c :::		   ghost region).
c ::: -----------------------------------------------------------
      subroutine FORT_INITDATA(level,time,lo,hi,nscal,
     &	 	                 vel,scal,DIMS(state),press,DIMS(press),
     &                         delta,xlo,xhi)
      implicit none
      integer    level, nscal, nvel
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     time, delta(SDIM)
      REAL_T     vel(DIMV(state),SDIM)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))

#include "htdata.H"
#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "INFL_FORCE_F.H"

c     ::::: local variables
      integer i, j, k, n, iN2, iO2
      REAL_T x, y, z, zt, u, v, w, rho, Yl(maxspec), Xl(maxspec), T, h, Patm
      REAL_T turbSclX, turbSclY, turbSclZ, convVelSAVE
      REAL_T eta, eta0, pmf_vals(maxspec+3), bg_vals(maxspec+3), r, y1, y2, theta
      character*(maxspnml) name

#if defined(BL_DO_FLCT)
      integer DIMDEC(uflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T uflct(:,:,:)
      allocatable uflct
      REAL_T fileDataZLO, fileDataZHI, gridLo(SDIM)
      integer bc(SDIM,2)

      do i = 1, SDIM
         gridLo(i) = xlo(i) - half*probSizeFile(i)
      end do
      fileDataZLO = half*(domnhi(3)+domnlo(3)-probSizeFile(3))
      fileDataZHI = half*(domnhi(3)+domnlo(3)+probSizeFile(3))
         
c     Make a platter of data in XY, will march over z
      do i = 1, SDIM
         loFlctArray(i) = lo(i)
         hiFlctArray(i) = hi(i)
      enddo

      strmwse_dir = 3
      loFlctArray(strmwse_dir) = 1
      hiFlctArray(strmwse_dir) = 1
      call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
      allocate(uflct(DIMV(uflct)))

      do i = 1, SDIM
         bc(i,1) = EXT_DIR
         bc(i,2) = EXT_DIR
      end do
#endif 
           
      do n = 1,Nspec
         call get_spec_name(name, n)
         if (name .eq. 'N2' ) iN2 = n
         if (name .eq. 'O2' ) iO2 = n
      end do

      convVelSAVE = convVel
      convVel = 15.d0
      do k= lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do n=1,SDIM
                  vel(i,j,k,n) = zero
               end do
            end do
         end do
            
         z = (float(k)+.5)*delta(3)+domnlo(3)
#if defined(BL_DO_FLCT)
c         if (.false. .and. flct_file.ne."") then      
         if (flct_file.ne."") then      
            call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo,
     &           delta, z, bc, domnlo, domnhi)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
#ifdef SWIRL
                  r = SQRT( ((float(i)+.5)*delta(1)+domnlo(1))**2
     &                 +    ((float(j)+.5)*delta(2)+domnlo(2))**2)
                  vel(i,j,k,1) = uflct(i,j,1)
     &                 * half*(one - TANH(two*(r-.015d0)/.002d0))
#else
                  vel(i,j,k,1) = uflct(i,j,1)*turb_scale
#endif
               end do
            end do
                  
            call INFL_FILL(FLCT_YVEL, DIMS(uflct), uflct, xlo,
     &           delta, z, bc, domnlo, domnhi)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
#ifdef SWIRL
                  r = SQRT( ((float(i)+.5)*delta(1)+domnlo(1))**2
     &                 +    ((float(j)+.5)*delta(2)+domnlo(2))**2)
                  vel(i,j,k,2) = uflct(i,j,1)
     &                 * half*(one - TANH(two*(r-.015d0)/.002d0))
#else
                  vel(i,j,k,2) = uflct(i,j,1)*turb_scale
#endif
               end do
            end do
               
            call INFL_FILL(FLCT_ZVEL, DIMS(uflct), uflct, xlo,
     &           delta, z, bc, domnlo, domnhi)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
#ifdef SWIRL
                  r = SQRT( ((float(i)+.5)*delta(1)+domnlo(1))**2
     &                 +    ((float(j)+.5)*delta(2)+domnlo(2))**2)
                  vel(i,j,k,3) = uflct(i,j,1)
     &                 * half*(one - TANH(two*(r-.015d0)/.002d0))
#else
                  vel(i,j,k,3) = uflct(i,j,1)*turb_scale
#endif
               end do
            end do
         end if
#endif
      enddo            
      convVel = convVelSAVE

      if (probtype .eq. 1) then
         do k = lo(3), hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
c
c Get V values from turb and lower bc, set Y and T separately
c
                  zt = domnlo(3) - half*delta(3)
                  call bcfunction(XLO,x,y,zt,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  if (z-domnlo(3) .GT. zstandoff) then
                     if ((SQRT(x*x+y*y).le.rhot).and.(z-domnlo(3).le.zhot+zstandoff)) T = Thot
                     do n = 1,Nspec
                        Xl(n) = zero
                     end do
                     Xl(iN2) = 0.79d0
                     Xl(iO2) = 0.21d0
                     CALL CKXTY (Xl, IWRK, RWRK, Yl)
                  end if
                  do n = 1,Nspec
                     scal(i,j,k,FirstSpec+n-1) = Yl(n)
                  end do
                  scal(i,j,k,Trac) = zero
                  scal(i,j,k,Temp) = T
#if 0
                  if (z-domnlo(3) .LT. zhot+zstandoff) then
                      vel(i,j,k,1) = u + vel(i,j,k,1)*turbSclX(x,y,z)
                      vel(i,j,k,2) = v + vel(i,j,k,2)*turbSclY(x,y,z)
                      vel(i,j,k,3) = w + vel(i,j,k,3)*turbSclZ(x,y,z)
                  else
                     vel(i,j,k,1) = zero
                     vel(i,j,k,2) = zero
                     vel(i,j,k,3) = zero
                  end if
#endif
               end do
            end do
         end do
      else if (probtype .eq. 2) then
         bg_vals(1) = T_bc(BL_VOLUME)
         do n = 1,Nspec
            bg_vals(3+n) = zero
         end do
         bg_vals(3+iN2) = 0.79d0
         bg_vals(3+iO2) = 0.21d0
         do k= lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)            
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  r = SQRT(x*x + y*y)

                  y1 = (z-half*delta(3)-zstandoff)*100.
                  y2 = (z+half*delta(3)-zstandoff)*100.
                  if (Nspec.eq.2) then
                     do n = 1,Nspec+3
                        pmf_vals(n) = bg_vals(n)
                     enddo
                  else
                     call pmf(y1,y2,pmf_vals,n)
                     if (n.ne.Nspec+3) then
                        write(6,*) 'INITDATA:  Wrong no. of species'
                        write(6,*) '           pmf file has  ',n-3
                        write(6,*) '           chem file has ',Nspec
                     endif
                  endif

                  eta = half*(one + TANH(two*(r-rhot)/zBL))

                  scal(i,j,k,Temp) = (one - eta)*pmf_vals(1) + eta*bg_vals(1)
                  do n = 1,Nspec
                     Xl(n) = (one - eta)*pmf_vals(3+n) + eta*bg_vals(3+n)
                  end do 

                  eta0 = half*(one + TANH(two*(z-zhot)/zBL))

                  scal(i,j,k,Temp) = (one - eta0)*scal(i,j,k,Temp) + eta0*bg_vals(1)
                  do n = 1,Nspec
                     Xl(n) = (one - eta0)*Xl(n) + eta0*bg_vals(3+n)
                  end do 

                  CALL CKXTY (Xl, IWRK, RWRK, Yl)
                  do n = 1,Nspec
                     scal(i,j,k,FirstSpec+n-1) = Yl(n)
                  end do
                  scal(i,j,k,Trac) = zero

                  zt = domnlo(3) - half*delta(3)

                  if (z-domnlo(3) .LT. vheight) then
                     call bcfunction(XLO,x,y,zt,time,u,v,w,rho,Yl,T,h,delta,.true.)
                      vel(i,j,k,1) = u + vel(i,j,k,1)*turbSclX(x,y,z)
                      vel(i,j,k,2) = v + vel(i,j,k,2)*turbSclY(x,y,z)
                      vel(i,j,k,3) = w + vel(i,j,k,3)*turbSclZ(x,y,z)
                  else
                     vel(i,j,k,1) = zero
                     vel(i,j,k,2) = zero
                     vel(i,j,k,3) = zero
                  end if
               end do
            end do
         end do
      else if (probtype .eq. 3) then
         bg_vals(1) = T_bc(BL_VOLUME)
         do n = 1,Nspec
            bg_vals(3+n) = zero
         end do
         bg_vals(3+iN2) = 0.79d0
         bg_vals(3+iO2) = 0.21d0
         theta = Pi / 6.d0
         do k= lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)            
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  r = SQRT(x*x + y*y)

                  if (stTh.ge.0.d0) then
                     y1 = ((z-ABS(x)*TAN(theta))-half*delta(3)-zstandoff)*100.
                  else
                     y1 = (z-zstandoff)*100.
                  end if

                  if (Nspec.eq.2) then
                     do n = 1,Nspec+3
                        pmf_vals(n) = bg_vals(n)
                     enddo
                  else
                     call pmf(y1,y1,pmf_vals,n)
                     if (n.ne.Nspec+3) then
                        write(6,*) 'INITDATA:  Wrong no. of species'
                        write(6,*) '           pmf file has  ',n-3
                        write(6,*) '           chem file has ',Nspec
                     endif
                  endif

                  eta = half*(one + TANH(two*(r-rhot)/zBL))

                  scal(i,j,k,Temp) = (one - eta)*pmf_vals(1) + eta*bg_vals(1)
                  do n = 1,Nspec
                     Xl(n) = (one - eta)*pmf_vals(3+n) + eta*bg_vals(3+n)
                  end do 

                  CALL CKXTY (Xl, IWRK, RWRK, Yl)
                  do n = 1,Nspec
                     scal(i,j,k,FirstSpec+n-1) = Yl(n)
                  end do
                  scal(i,j,k,Trac) = zero

                  zt = domnlo(3) - half*delta(3)
                  if (z-domnlo(3) .LT. vheight) then
                     vel(i,j,k,1) = zero
                     vel(i,j,k,2) = zero
                     vel(i,j,k,3) = zero
                  end if
               end do
            end do
         end do
      end if

      Patm = pamb / 101325.0d0
      call FORT_RHOfromPTY(lo,hi,
     &     scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Density),  DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state),
     &     Patm)
      call FORT_HMIXfromTY(lo,hi,
     &     scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),RhoH),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state))

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do n = 0,Nspec-1
                  scal(i,j,k,FirstSpec+n) = scal(i,j,k,FirstSpec+n)*scal(i,j,k,Density)
               enddo
               scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
            enddo
         enddo
      enddo

#if defined(BL_DO_FLCT)          
      deallocate(uflct)
#endif
      end
      
c ::: -----------------------------------------------------------
c ::: This routine will zero out diffusivity on portions of the
c ::: boundary that are inflow, allowing that a "wall" block
c ::: the complement aperture
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: diff      <=> diffusivity on edges
c ::: DIMS(diff) => index extent of diff array
c ::: lo,hi      => region of interest, edge-based
c ::: domlo,hi   => index extent of problem domain, edge-based
c ::: dx         => cell spacing
c ::: problo     => phys loc of lower left corner of prob domain
c ::: bc         => boundary condition flag (on orient)
c :::                   in BC_TYPES::physicalBndryTypes
c ::: idir       => which face, 0=x, 1=y
c ::: isrz       => 1 if problem is r-z
c ::: id         => index of state, 0=u
c ::: ncomp      => components to modify
c ::: 
c ::: -----------------------------------------------------------
      subroutine FORT_ZEROVISC(diff,DIMS(diff),lo,hi,domlo,domhi,
     &                         dx,problo,bc,idir,isrz,id,ncomp)
      implicit none
      integer DIMDEC(diff)
      integer lo(SDIM), hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      integer bc(2*SDIM)
      integer idir, isrz, id, ncomp
      REAL_T  diff(DIMV(diff),*)
      REAL_T  dx(SDIM)
      REAL_T  problo(SDIM)
      
#include "probdata.H"
#include "htdata.H"
      integer i, j, k, n, Tid, RHid, YSid, YEid, ys, ye
      integer getZone
      logical do_T, do_RH, do_Y
      REAL_T xl, xr, xh, yb, yt, yh, z

      Tid  = Temp      - id + SDIM
      RHid = RhoH      - id + SDIM
      YSid = FirstSpec - id + SDIM
      YEid = LastSpec  - id + SDIM
      
      do_T  = (Tid  .GE. 1) .AND. (Tid  .LE. ncomp)
      do_RH = (RHid .GE. 1) .AND. (RHid .LE. ncomp)
      ys = MAX(YSid,1)
      ye = MIN(YEid,ncomp)
      do_Y = (ye - ys + 1) .GE. 1

c     Do species, Temp, rhoH
      if ((idir.EQ.2) .AND. (lo(3) .LE. domlo(3))
     &        .AND. (do_T .OR. do_RH .OR. do_Y) ) then

         k = lo(3)
         z = float(k)*dx(3)+domnlo(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               xl = float(i)*dx(1)+domnlo(1) 
               xr = (float(i)+1.d0)*dx(1)+domnlo(1) 
               xh = 0.5d0*(xl+xh)
               yb = float(j)*dx(2)+domnlo(2) 
               yt = (float(j)+1.d0)*dx(2)+domnlo(2) 
               yh = 0.5d0*(yb+yt)
               
               if ( (getZone(xl,yb,z).eq.BL_STICK) .OR.
     &              (getZone(xh,yb,z).eq.BL_STICK) .OR.
     &              (getZone(xr,yb,z).eq.BL_STICK) .OR.
     &              (getZone(xl,yh,z).eq.BL_STICK) .OR.
     &              (getZone(xh,yh,z).eq.BL_STICK) .OR.
     &              (getZone(xr,yh,z).eq.BL_STICK) .OR.
     &              (getZone(xl,yt,z).eq.BL_STICK) .OR.
     &              (getZone(xh,yt,z).eq.BL_STICK) .OR.
     &              (getZone(xr,yt,z).eq.BL_STICK) .OR.
     &              (getZone(xl,yb,z).eq.BL_PIPEEND) .OR.
     &              (getZone(xh,yb,z).eq.BL_PIPEEND) .OR.
     &              (getZone(xr,yb,z).eq.BL_PIPEEND) .OR.
     &              (getZone(xl,yh,z).eq.BL_PIPEEND) .OR.
     &              (getZone(xh,yh,z).eq.BL_PIPEEND) .OR.
     &              (getZone(xr,yh,z).eq.BL_PIPEEND) .OR.
     &              (getZone(xl,yt,z).eq.BL_PIPEEND) .OR.
     &              (getZone(xh,yt,z).eq.BL_PIPEEND) .OR.
     &              (getZone(xr,yt,z).eq.BL_PIPEEND) ) then

                  if (do_T)  diff(i,j,k,Tid ) = 0.d0
                  if (do_RH) diff(i,j,k,RHid) = 0.d0
                  if (do_Y) then
                     do n=ys,ye
                        diff(i,j,k,n) = 0.d0
                     enddo
                  endif

               endif
            end do
         end do
      end if
      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: density gradient
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: rho       => density array
c ::: DIMS(rho) => index extent of rho array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_DENERROR (tag,DIMS(tag),set,clear,
     &                          rho,DIMS(rho),lo,hi,nvar,
     &                          domlo,domhi,dx,xlo,
     &			        problo,time,level)
      implicit none
      integer   DIMDEC(rho)
      integer   DIMDEC(tag)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    rho(DIMV(rho), nvar)

      integer   i, j

#include "probdata.H"

      call bl_abort('DENERROR: should no be here for H-flame problem')
      
      end

c ::: -----------------------------------------------------------

      subroutine FORT_FLAMETRACERROR (tag,DIMS(tag),set,clear,
     &                                ftrac,DIMS(ftrac),lo,hi,nvar,
     &                                domlo,domhi,dx,xlo,
     &			              problo,time,level)
      implicit none
      integer   DIMDEC(ftrac)
      integer   DIMDEC(tag)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    ftrac(DIMV(ftrac), nvar)

      integer   i, j, k
      REAL_T    x, y, z, r
      logical   in_refine_zone

#include "probdata.H"

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = cvmgt(set,tag(i,j,k),
     &                            ftrac(i,j,k,1).gt.flametracval)
            enddo
         enddo
      enddo

      if ((level.lt.max_nozzle_lev).and.(refine_nozzle.eq.1)) then
         do k = lo(3), hi(3)
            z = (float(k)+.5)*dx(3)+problo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*dx(2)+problo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*dx(1)+problo(1)
                  r = SQRT(x*x+y*y)
                  in_refine_zone = (r .le. refine_nozzle_x)
                  in_refine_zone = in_refine_zone .and.
     &                  (z - problo(3)) .le. refine_nozzle_z
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),in_refine_zone)
               end do
            end do
         end do
      end if
      if ((level.lt.max_stick_lev).and.(refine_stick.eq.1)) then
         do k = lo(3), hi(3)
            z = (float(k)+.5)*dx(3)+problo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*dx(2)+problo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*dx(1)+problo(1)
                  r = SQRT(x*x+y*y)
                  in_refine_zone = ((ABS(x).le.refine_stick_x).and.(r.le.Ro))
                  in_refine_zone = in_refine_zone .and.
     &                  (z - problo(3)) .le. refine_stick_z
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),in_refine_zone)
               end do
            end do
         end do
      end if
      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: density gradient
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: adv       => scalar array
c ::: DIMS(adv) => index extent of scalar array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_ADVERROR (tag,DIMS(tag),set,clear,
     &                          adv,DIMS(adv),lo,hi,nvar,
     &                          domlo,domhi,delta,xlo,
     &			        problo,time,level)
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      REAL_T    delta(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    adv(DIMV(adv),nvar)

#include "probdata.H"
      
      integer i, j

      if ((probtype .eq. 8 .or. probtype .eq. 28 .or. probtype .eq. 29) 
     &     .and. (time .eq. zero)) then
         call FORT_MVERROR(tag,DIMS(tag),set,clear,
     &                     adv,DIMS(adv),lo,hi,nvar,
     &                     domlo,domhi,delta,xlo,
     &                     problo,time,level)
      endif
      
      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the
c ::: temperature gradient
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: temp      => density array
c ::: DIMS(temp)=> index extent of temp array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::              corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_TEMPERROR (tag,DIMS(tag),set,clear,
     &                          temperature,DIMS(temp),lo,hi,nvar,
     &                          domlo,domhi,dx,xlo,
     &                          problo,time,level)
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(temp)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    temperature(DIMV(temp),nvar)

      REAL_T    x, y, z, ax, ay, az, aerr
      integer   i, j, k, ng

#include "probdata.H"

      ng = min(ARG_H1(temp)-hi(1),ARG_H2(temp)-hi(2),ARG_H3(temp)-hi(3),
     &         lo(1)-ARG_L1(temp),lo(2)-ARG_L2(temp),lo(3)-ARG_L3(temp))

      if (ng .lt. 1) then
         write(6,*) "TEMPERR cannot compute gradient, ng = ",ng
         call bl_abort(" ")
      endif

c     ::::: refine where there is temperature gradient
      if (level .lt. max_temp_lev) then
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ax = abs(temperature(i+1,j,k,1) - temperature(i,j,k,1))
            ay = abs(temperature(i,j+1,k,1) - temperature(i,j,k,1))
            az = abs(temperature(i,j,k+1,1) - temperature(i,j,k,1))
            ax = MAX(ax,abs(temperature(i,j,k,1) - temperature(i-1,j,k,1)))
            ay = MAX(ay,abs(temperature(i,j,k,1) - temperature(i,j-1,k,1)))
            az = MAX(az,abs(temperature(i,j,k,1) - temperature(i,j,k-1,1)))
            aerr = max(ax,ay,az)
            tag(i,j,k) = cvmgt(set,tag(i,j,k),aerr.ge.tempgrad)
         enddo
      enddo
      enddo
      endif

      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: magnitude of vorticity
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: vort      => array of vorticity values
c ::: DIMS(vor) => index extent of vort array
c ::: nvar      => number of components in vort array (should be 1)
c ::: lo,hi     => index extent of grid
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_MVERROR (tag,DIMS(tag),set,clear,
     &                         vort,DIMS(vort),lo,hi,nvar,
     &                         domlo,domhi,dx,xlo,
     &			       problo,time,level)
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(vort)
      integer   nvar, set, clear, level
      integer   lo(SDIM), hi(SDIM)
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    vort(DIMV(vort),nvar)

      REAL_T    x, y ,z
      integer   i, j, k
      REAL_T    max_vort, tagval

#include "probdata.H"

c     if ((probtype.eq.8 .or. probtype .eq. 28 .or. probtype .eq. 29) 
c    &     .and. (level .lt. max_vort_lev)) then
      if (level .lt. max_vort_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),
     &                         ABS(vort(i,j,k,1)).ge.vorterr)
               enddo
            enddo
         enddo
      end if

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: den      <=  density array
c ::: DIMS(den) => index extent of den array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of den array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_DENFILL (den,DIMS(den),domlo,domhi,delta,
     &                         xlo,time,bc)
      implicit none

      integer DIMDEC(den), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  den(DIMV(den))

#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(den)
      lo(2) = ARG_L2(den)
      lo(3) = ARG_L3(den)
      hi(1) = ARG_H1(den)
      hi(2) = ARG_H2(den)
      hi(3) = ARG_H3(den)

      ilo = max(lo(1),domlo(1))
      jlo = max(lo(2),domlo(2))
      klo = max(lo(3),domlo(3))
      ihi = min(hi(1),domhi(1))
      jhi = min(hi(2),domhi(2))
      khi = min(hi(3),domhi(3))
      
      call filcc (den,DIMS(den),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif    
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: adv      <=  advected quantity array
c ::: DIMS(adv) => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of adv array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,delta,xlo,time,bc)

      implicit none

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)

      call filcc (adv,DIMS(adv),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(adv,DIMS(adv),domlo,domhi,bc)

      end


c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.
c :::
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: temp     <=  temperature array
c ::: lo,hi     => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::              corner of temperature array
c ::: time      => problem evolution time
c ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_TEMPFILL (temp,DIMS(temp),domlo,domhi,delta,
     &                         xlo,time,bc)

      implicit none

      integer DIMDEC(temp), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  temp(DIMV(temp))

#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(temp)
      lo(2) = ARG_L2(temp)
      lo(3) = ARG_L3(temp)
      hi(1) = ARG_H1(temp)
      hi(2) = ARG_H2(temp)
      hi(3) = ARG_H3(temp)

      ilo = max(lo(1),domlo(1))
      jlo = max(lo(2),domlo(2))
      klo = max(lo(3),domlo(3))
      ihi = min(hi(1),domhi(1))
      jhi = min(hi(2),domhi(2))
      khi = min(hi(3),domhi(3))
      
      call filcc (temp,DIMS(temp),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif    
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.
c :::
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: rhoh      <=  rho*h array
c ::: lo,hi     => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::              corner of temperature array
c ::: time      => problem evolution time
c ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_RHOHFILL (rhoh,DIMS(rhoh),domlo,domhi,delta,
     &                         xlo,time,bc)

      implicit none

      integer DIMDEC(rhoh), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoh(DIMV(rhoh))

#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(rhoh)
      lo(2) = ARG_L2(rhoh)
      lo(3) = ARG_L3(rhoh)
      hi(1) = ARG_H1(rhoh)
      hi(2) = ARG_H2(rhoh)
      hi(3) = ARG_H3(rhoh)

      ilo = max(lo(1),domlo(1))
      jlo = max(lo(2),domlo(2))
      klo = max(lo(3),domlo(3))
      ihi = min(hi(1),domhi(1))
      jhi = min(hi(2),domhi(2))
      khi = min(hi(3),domhi(3))
      
      call filcc (rhoh,DIMS(rhoh),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif    
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif

      end
c
c Fill x, y & z velocity at once.
c
      subroutine FORT_VELFILL (vel,DIMS(vel),domlo,domhi,delta,
     &                         xlo,time,bc)

      implicit none

      integer DIMDEC(vel), bc(SDIM,2,SDIM)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  vel(DIMV(vel),SDIM)

#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "INFL_FORCE_F.H"

      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
      REAL_T  cs,sg,r,scaler,scalet
      REAL_T  turbSclX,turbSclY,turbSclZ,turbSclR,turbSclT

#if defined(BL_DO_FLCT)
      integer DIMDEC(uflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T  uflct(:,:,:), vflct(:,:,:),  wflct(:,:,:)
      allocatable uflct, vflct, wflct
#endif
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(vel)
      lo(2) = ARG_L2(vel)
      lo(3) = ARG_L3(vel)
      hi(1) = ARG_H1(vel)
      hi(2) = ARG_H2(vel)
      hi(3) = ARG_H3(vel)

      ilo = max(lo(1),domlo(1))
      jlo = max(lo(2),domlo(2))
      klo = max(lo(3),domlo(3))
      ihi = min(hi(1),domhi(1))
      jhi = min(hi(2),domhi(2))
      khi = min(hi(3),domhi(3)) 

#if defined(BL_DO_FLCT)
      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         enddo
         loFlctArray(strmwse_dir) = 1
         hiFlctArray(strmwse_dir) = 1
         call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
         allocate(uflct(DIMV(uflct)))
         allocate(vflct(DIMV(uflct)))
         allocate(wflct(DIMV(uflct)))
         call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo, delta, time,
     $                  bc(1,1,1), domnlo, domnhi)
         call INFL_FILL(FLCT_YVEL, DIMS(uflct), vflct, xlo, delta, time,
     $                  bc(1,1,2), domnlo, domnhi)
         call INFL_FILL(FLCT_ZVEL, DIMS(uflct), wflct, xlo, delta, time,
     $                  bc(1,1,3), domnlo, domnhi)
      endif
#endif

      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),1),
     $            DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,1))
      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),2),
     $            DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,2))
      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),3),
     $            DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,3))

      do i = lo(1), domlo(1)-1
         x = (float(i)+.5)*delta(1)+domnlo(1)
         do k = lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)

               call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

               if (bc(1,1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 1) then
                     vel(i,j,k,1) = u + uflct(1,j,k)*turbSclX(x,y,z)*turb_scale
                  else
                     vel(i,j,k,1) = u
                  endif
#else
                  vel(i,j,k,1) = u
#endif
               endif

               if (bc(1,1,2).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 1) then
                     vel(i,j,k,2) = v + vflct(1,j,k)*turbSclY(x,y,z)*turb_scale
                  else
                     vel(i,j,k,2) = v
                  endif
#else
                  vel(i,j,k,2) = v
#endif
               endif

               if (bc(1,1,3).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 1) then
                     vel(i,j,k,3) = w + wflct(1,j,k)*turbSclZ(x,y,z)*turb_scale
                  else
                     vel(i,j,k,3) = w
                  endif
#else
                  vel(i,j,k,3) = w
#endif
               endif

            enddo
         enddo
      enddo

      do i = domhi(1)+1, hi(1)
         x = (float(i)+.5)*delta(1)+domnlo(1)
         do k = lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)

               call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

               if (bc(1,2,1).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 1) then
                     vel(i,j,k,1) = u + uflct(1,j,k)*turbSclX(x,y,z)*turb_scale
                  else
                     vel(i,j,k,1) = u
                  endif
#else
                  vel(i,j,k,1) = u
#endif
               endif

               if (bc(1,2,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 1) then
                     vel(i,j,k,2) = v + vflct(1,j,k)*turbSclY(x,y,z)*turb_scale
                  else
                     vel(i,j,k,2) = v
                  endif
#else
                  vel(i,j,k,2) = v
#endif
               endif

               if (bc(1,2,3).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 1) then
                     vel(i,j,k,3) = w + wflct(1,j,k)*turbSclZ(x,y,z)*turb_scale
                  else
                     vel(i,j,k,3) = w
                  endif
#else
                  vel(i,j,k,3) = w
#endif
               endif

            enddo
         enddo
      enddo

      do j = lo(2), domlo(2)-1
         y = (float(j)+.5)*delta(2)+domnlo(2)
         do k = lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)

               call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

               if (bc(2,1,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 2) then
                     vel(i,j,k,1) = u + uflct(i,1,k)*turbSclX(x,y,z)*turb_scale
                  else
                     vel(i,j,k,1) = u
                  endif
#else
                  vel(i,j,k,1) = u
#endif
               endif

               if (bc(2,1,2).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 2) then
                     vel(i,j,k,2) = v + vflct(i,1,k)*turbSclY(x,y,z)*turb_scale
                  else
                     vel(i,j,k,2) = v
                  endif
#else
                  vel(i,j,k,2) = v
#endif
               endif

               if (bc(2,1,3).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 2) then
                     vel(i,j,k,3) = w + wflct(i,1,k)*turbSclZ(x,y,z)*turb_scale
                  else
                     vel(i,j,k,3) = w
                  endif
#else
                  vel(i,j,k,3) = w
#endif
               endif

            enddo
         enddo
      enddo

      do j = domhi(2)+1, hi(2)
         y = (float(j)+.5)*delta(2)+domnlo(2)
         do k = lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)

               call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

               if (bc(2,2,1).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 2) then
                     vel(i,j,k,1) = u + uflct(i,1,k)*turbSclX(x,y,z)*turb_scale
                  else
                     vel(i,j,k,1) = u
                  endif
#else
                  vel(i,j,k,1) = u
#endif
               endif

               if (bc(2,2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 2) then
                     vel(i,j,k,2) = v + vflct(i,1,k)*turbSclY(x,y,z)*turb_scale
                  else
                     vel(i,j,k,2) = v
                  endif
#else
                  vel(i,j,k,2) = v
#endif
               endif

               if (bc(2,2,3).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 2) then
                     vel(i,j,k,3) = w + wflct(i,1,k)*turbSclZ(x,y,z)*turb_scale
                  else
                     vel(i,j,k,3) = w
                  endif
#else
                  vel(i,j,k,3) = w
#endif
               endif

            enddo
         enddo
      enddo

      do k = lo(3), domlo(3)-1
         z = (float(k)+.5)*delta(3)+domnlo(3)
         do j = lo(2),hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)

               call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

               if (bc(3,1,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 3) then
#ifdef SWIRL
                     vel(i,j,k,1) = u + uflct(i,j,1)
#else
                     r = sqrt(x**2+y**2)
                     cs = x/r
                     sg = y/r
                     scaler = turbSclR(r)
                     scalet = turbSclT(r)
                     vel(i,j,k,1) = u + turb_scale*turbSclX(x,y,z)*
     &                 ((scaler*cs**2+scalet*sg**2)*uflct(i,j,1) +
     &                  (scaler*sg*cs-scalet*sg*cs)*vflct(i,j,1))
#endif
                  else
                     vel(i,j,k,1) = u
                  endif
#else
                  vel(i,j,k,1) = u
#endif
               endif

               if (bc(3,1,2).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 3) then
#ifdef SWIRL
                     vel(i,j,k,2) = v + vflct(i,j,1)
#else
                     r = sqrt(x**2+y**2)
                     cs = x/r
                     sg = y/r
                     scaler = turbSclR(r)
                     scalet = turbSclT(r)
                     vel(i,j,k,2) = v + turb_scale*turbSclY(x,y,z)*
     &                 ((scaler*sg**2+scalet*cs**2)*vflct(i,j,1) +
     &                  (scaler*sg*cs- scalet*sg*cs)*uflct(i,j,1))
#endif
                  else
                     vel(i,j,k,2) = v
                  endif
#else
                  vel(i,j,k,2) = v
#endif
               endif

               if (bc(3,1,3).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 3) then
#ifdef SWIRL
                     vel(i,j,k,3) = w + wflct(i,j,1)
#else
                     vel(i,j,k,3) = w + wflct(i,j,1)*turbSclZ(x,y,z)*turb_scale
     &                      * anisotsc
#endif
                  else
                     vel(i,j,k,3) = w
                  endif
#else
                  vel(i,j,k,3) = w
#endif
               endif

            enddo
         enddo
      enddo

      do k = domhi(3)+1, hi(3)
         z = (float(k)+.5)*delta(3)+domnlo(3)
         do j = lo(2),hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)

               call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

               if (bc(3,2,1).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 3) then
                     vel(i,j,k,1) = u + uflct(i,j,1)*turbSclX(x,y,z)*turb_scale
                  else
                     vel(i,j,k,1) = u
                  endif
#else
                  vel(i,j,k,1) = u
#endif
               endif

               if (bc(3,2,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 3) then
                     vel(i,j,k,2) = v + vflct(i,j,1)*turbSclY(x,y,z)*turb_scale
                  else
                     vel(i,j,k,2) = v
                  endif
#else
                  vel(i,j,k,2) = v
#endif
               endif

               if (bc(3,2,3).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 3) then
                     vel(i,j,k,3) = w + wflct(i,j,1)*turbSclZ(x,y,z)*turb_scale
                  else
                     vel(i,j,k,3) = w
                  endif
#else
                  vel(i,j,k,3) = w
#endif
               endif

            enddo
         enddo
      enddo


#if defined(BL_DO_FLCT)
      if (forceInflow) then
         deallocate(uflct)
         deallocate(vflct)
         deallocate(wflct)
      endif
#endif

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: xvel     <=  x velocity array
c ::: lo,hi     => index extent of xvel array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_XVELFILL (xvel,DIMS(xvel),domlo,domhi,delta,
     &                          xlo,time,bc)

      implicit none

      integer DIMDEC(xvel), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  xvel(DIMV(xvel))

#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "INFL_FORCE_F.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
      REAL_T  cs,sg,r,scaler,scalet
      REAL_T turbSclX,turbSclR,turbSclT
#if defined(BL_DO_FLCT)
      integer DIMDEC(uflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T uflct(:,:,:), vflct(:,:,:)
      allocatable uflct, vflct
#endif
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(xvel)
      lo(2) = ARG_L2(xvel)
      lo(3) = ARG_L3(xvel)
      hi(1) = ARG_H1(xvel)
      hi(2) = ARG_H2(xvel)
      hi(3) = ARG_H3(xvel)

      ilo = max(lo(1),domlo(1))
      jlo = max(lo(2),domlo(2))
      klo = max(lo(3),domlo(3))
      ihi = min(hi(1),domhi(1))
      jhi = min(hi(2),domhi(2))
      khi = min(hi(3),domhi(3)) 

#if defined(BL_DO_FLCT)
      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         enddo
         loFlctArray(strmwse_dir) = 1
         hiFlctArray(strmwse_dir) = 1
         call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
         allocate(uflct(DIMV(uflct)))
         call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo, delta, time,
     $                  bc, domnlo, domnhi)
         allocate(vflct(DIMV(uflct)))
         call INFL_FILL(FLCT_YVEL, DIMS(uflct), vflct, xlo, delta, time,
     $                  bc, domnlo, domnhi)
      endif
#endif

      call filcc (xvel,DIMS(xvel),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 1) then
                     xvel(i,j,k) = u + uflct(1,j,k)*turbSclX(x,y,z)*turb_scale
                  else
                     xvel(i,j,k) = u
                  endif
#else
                  xvel(i,j,k) = u
#endif
               enddo
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 1) then
                     xvel(i,j,k) = u + uflct(1,j,k)*turbSclX(x,y,z)*turb_scale
                  else
                     xvel(i,j,k) = u
                  endif
#else
                  xvel(i,j,k) = u
#endif
               enddo
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 2) then
                     xvel(i,j,k) = u + uflct(i,1,k)*turbSclX(x,y,z)*turb_scale
                  else
                     xvel(i,j,k) = u
                  endif
#else
                  xvel(i,j,k) = u
#endif
               enddo
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 2) then
                     xvel(i,j,k) = u + uflct(i,1,k)*turbSclX(x,y,z)*turb_scale
                  else
                     xvel(i,j,k) = u
                  endif
#else
                  xvel(i,j,k) = u
#endif
               enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 3) then
                     r = sqrt(x**2+y**2)
                     cs = x/r
                     sg = y/r
                     scaler = turbSclR(r)
                     scalet = turbSclT(r)
                     xvel(i,j,k) = u + turb_scale*turbSclX(x,y,z)*
     &                 ((scaler*cs**2+scalet*sg**2)*uflct(i,j,1) +
     &                  (scaler*sg*cs-scalet*sg*cs)*vflct(i,j,1))
                  else
                     xvel(i,j,k) = u
                  endif
#else
                  xvel(i,j,k) = u
#endif
               enddo
            enddo
         enddo
      endif    
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 3) then
                     xvel(i,j,k) = u + uflct(i,j,1)*turbSclX(x,y,z)*turb_scale
                  else
                     xvel(i,j,k) = u
                  endif
#else
                  xvel(i,j,k) = u
#endif
               enddo
            enddo
         enddo
      endif

#if defined(BL_DO_FLCT)
      if (forceInflow) then
         deallocate(uflct)
         deallocate(vflct)
      endif
#endif

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: yvel     <=  y velocity array
c ::: lo,hi     => index extent of yvel array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_YVELFILL (yvel,DIMS(yvel),domlo,domhi,delta,
     &                          xlo,time,bc)

      implicit none

      integer DIMDEC(yvel), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  yvel(DIMV(yvel))

#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "INFL_FORCE_F.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
      REAL_T  cs,sg,r,scaler,scalet
      REAL_T turbSclY,turbSclR,turbSclT

#if defined(BL_DO_FLCT)
      integer DIMDEC(vflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T uflct(:,:,:), vflct(:,:,:)
      allocatable uflct, vflct
#endif

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(yvel)
      lo(2) = ARG_L2(yvel)
      lo(3) = ARG_L3(yvel)
      hi(1) = ARG_H1(yvel)
      hi(2) = ARG_H2(yvel)
      hi(3) = ARG_H3(yvel)

      ilo = max(lo(1),domlo(1))
      jlo = max(lo(2),domlo(2))
      klo = max(lo(3),domlo(3))
      ihi = min(hi(1),domhi(1))
      jhi = min(hi(2),domhi(2))
      khi = min(hi(3),domhi(3))
      
#if defined(BL_DO_FLCT)
      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         enddo
         loFlctArray(strmwse_dir) = 1
         hiFlctArray(strmwse_dir) = 1
         call SET_ARGS(DIMS(vflct), loFlctArray, hiFlctArray)
         allocate(uflct(DIMV(vflct)))
         call INFL_FILL(FLCT_XVEL, DIMS(vflct), uflct, xlo, delta, time,
     $                  bc, domnlo, domnhi)
         allocate(vflct(DIMV(vflct)))
         call INFL_FILL(FLCT_YVEL, DIMS(vflct), vflct, xlo, delta, time,
     $                  bc, domnlo, domnhi)
      endif
#endif

      call filcc (yvel,DIMS(yvel),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 1) then
                     yvel(i,j,k) = v + vflct(1,j,k)*turbSclY(x,y,z)*turb_scale
                  else
                     yvel(i,j,k) = v
                  endif
#else
                  yvel(i,j,k) = v
#endif
               enddo
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 1) then
                     yvel(i,j,k) = v + vflct(1,j,k)*turbSclY(x,y,z)*turb_scale
                  else
                     yvel(i,j,k) = v
                  endif
#else
                  yvel(i,j,k) = v
#endif
               enddo
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 2) then
                     yvel(i,j,k) = v + vflct(i,1,k)*turbSclY(x,y,z)*turb_scale
                  else
                     yvel(i,j,k) = v
                  endif
#else
                  yvel(i,j,k) = v
#endif
               enddo
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 2) then
                     yvel(i,j,k) = v + vflct(i,1,k)*turbSclY(x,y,z)*turb_scale
                  else
                     yvel(i,j,k) = v
                  endif
#else
                  yvel(i,j,k) = v
#endif
               enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 3) then
c                    yvel(i,j,k) = v + vflct(i,j,1)*turbSclY(x,y,z)*turb_scale
                     r = sqrt(x**2+y**2)
                     cs = x/r
                     sg = y/r
                     scaler = turbSclR(r)
                     scalet = turbSclT(r)
                     yvel(i,j,k) = v + turb_scale*turbSclY(x,y,z)*
     &                 ((scaler*sg**2+scalet*cs**2)*vflct(i,j,1) +
     &                  (scaler*sg*cs- scalet*sg*cs)*uflct(i,j,1))
                  else
                     yvel(i,j,k) = v
                  endif
#else
                  yvel(i,j,k) = v
#endif
               enddo
            enddo
         enddo
      endif    
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 3) then
                     yvel(i,j,k) = v + vflct(i,j,1)*turbSclY(x,y,z)*turb_scale
                  else
                     yvel(i,j,k) = v
                  endif
#else
                  yvel(i,j,k) = v
#endif
               enddo
            enddo
         enddo
      endif

#if defined(BL_DO_FLCT)
      if (forceInflow) then
         deallocate(uflct)
         deallocate(vflct)
      endif
#endif

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: zvel     <=  z velocity array
c ::: lo,hi     => index extent of zvel array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ZVELFILL (zvel,DIMS(zvel),domlo,domhi,delta,
     &                          xlo,time,bc)

      implicit none

      integer DIMDEC(zvel), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  zvel(DIMV(zvel))

#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "INFL_FORCE_F.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
      REAL_T turbSclZ
#if defined(BL_DO_FLCT)
      integer DIMDEC(wflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T wflct(:,:,:)
      allocatable wflct
#endif

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(zvel)
      lo(2) = ARG_L2(zvel)
      lo(3) = ARG_L3(zvel)
      hi(1) = ARG_H1(zvel)
      hi(2) = ARG_H2(zvel)
      hi(3) = ARG_H3(zvel)

      ilo = max(lo(1),domlo(1))
      jlo = max(lo(2),domlo(2))
      klo = max(lo(3),domlo(3))
      ihi = min(hi(1),domhi(1))
      jhi = min(hi(2),domhi(2))
      khi = min(hi(3),domhi(3))

#if defined(BL_DO_FLCT)
      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         enddo
         loFlctArray(strmwse_dir) = 1
         hiFlctArray(strmwse_dir) = 1
         call SET_ARGS(DIMS(wflct), loFlctArray, hiFlctArray)
         allocate(wflct(DIMV(wflct)))
         call INFL_FILL(FLCT_ZVEL, DIMS(wflct), wflct, xlo, delta, time,
     $                  bc, domnlo, domnhi)
      endif
#endif
      
      call filcc (zvel,DIMS(zvel),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 1) then
                     zvel(i,j,k) = w + wflct(1,j,k)*turbSclZ(x,y,z)*turb_scale
                  else
                     zvel(i,j,k) = w
                  endif
#else
                  zvel(i,j,k) = w
#endif
               enddo
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 1) then
                     zvel(i,j,k) = w + wflct(1,j,k)*turbSclZ(x,y,z)*turb_scale
                  else
                     zvel(i,j,k) = w
                  endif
#else
                  zvel(i,j,k) = w
#endif
               enddo
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 2) then
                     zvel(i,j,k) = w + wflct(i,1,k)*turbSclZ(x,y,z)*turb_scale
                  else
                     zvel(i,j,k) = w
                  endif
#else
                  zvel(i,j,k) = w
#endif
               enddo
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 2) then
                     zvel(i,j,k) = w + wflct(i,1,k)*turbSclZ(x,y,z)*turb_scale
                  else
                     zvel(i,j,k) = w
                  endif
#else
                  zvel(i,j,k) = w
#endif
               enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceLo .and. strmwse_dir .eq. 3) then
                     zvel(i,j,k) = w + wflct(i,j,1)*turbSclZ(x,y,z)*turb_scale
     &                      * anisotsc
                  else
                     zvel(i,j,k) = w
                  endif
#else
                  zvel(i,j,k) = w
#endif
               enddo
            enddo
         enddo
      endif    
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                  if (forceHi .and. strmwse_dir .eq. 3) then
                     zvel(i,j,k) = w + wflct(i,j,1)*turbSclZ(x,y,z)*turb_scale
                  else
                     zvel(i,j,k) = w
                  endif
#else
                  zvel(i,j,k) = w
#endif
               enddo
            enddo
         enddo
      endif

#if defined(BL_DO_FLCT)
      if (forceInflow) then
         deallocate(wflct)
      endif
#endif

      end

      subroutine FORT_ALLCHEMFILL (rhoY,DIMS(rhoY),domlo,domhi,delta,
     &                             xlo,time,bc)

      implicit none

#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer DIMDEC(rhoY), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoY(DIMV(rhoY),Nspec)

      integer i, j, k, n
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(maxspec), T, h

      integer lo(SDIM), hi(SDIM)

c      print *, 'FORT_ALLCHEMFILL: ', domlo,domhi,delta,xlo,time

      lo(1) = ARG_L1(rhoY)
      lo(2) = ARG_L2(rhoY)
      lo(3) = ARG_L3(rhoY)
      hi(1) = ARG_H1(rhoY)
      hi(2) = ARG_H2(rhoY)
      hi(3) = ARG_H3(rhoY)

      ilo = max(lo(1),domlo(1))
      jlo = max(lo(2),domlo(2))
      klo = max(lo(3),domlo(3))
      ihi = min(hi(1),domhi(1))
      jhi = min(hi(2),domhi(2))
      khi = min(hi(3),domhi(3))
      
      do n = 1,Nspec
         call filcc (rhoY(lo(1),lo(2),lo(3),n),
     $               DIMS(rhoY),domlo,domhi,delta,xlo,bc)
      end do

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  do n = 1,Nspec
                     rhoY(i,j,k,n) = rho*Yl(n)
                  end do
               enddo
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  do n = 1,Nspec
                     rhoY(i,j,k,n) = rho*Yl(n)
                  end do
               enddo
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  do n = 1,Nspec
                     rhoY(i,j,k,n) = rho*Yl(n)
                  end do
               enddo
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  do n = 1,Nspec
                     rhoY(i,j,k,n) = rho*Yl(n)
                  end do
               enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  do n = 1,Nspec
                     rhoY(i,j,k,n) = rho*Yl(n)
                  end do
               enddo
            enddo
         enddo
      endif    
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  do n = 1,Nspec
                     rhoY(i,j,k,n) = rho*Yl(n)
                  end do
               enddo
            enddo
         enddo
      endif

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.
c :::
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: rhoY      <= rho*Y (Y=mass fraction) array
c ::: lo,hi     => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::              corner of temperature array
c ::: time      => problem evolution time
c ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: stateID   => id index of state being filled
c ::: -----------------------------------------------------------

      subroutine FORT_CHEMFILL (rhoY,DIMS(rhoY),domlo,domhi,delta,
     &                         xlo,time,bc,id)

      implicit none

      integer DIMDEC(rhoY), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM), id
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoY(DIMV(rhoY))

#include "ckdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

c      print *, 'FORT_CHEMFILL: ', domlo,domhi,delta,xlo,time

      lo(1) = ARG_L1(rhoY)
      lo(2) = ARG_L2(rhoY)
      lo(3) = ARG_L3(rhoY)
      hi(1) = ARG_H1(rhoY)
      hi(2) = ARG_H2(rhoY)
      hi(3) = ARG_H3(rhoY)

      ilo = max(lo(1),domlo(1))
      jlo = max(lo(2),domlo(2))
      klo = max(lo(3),domlo(3))
      ihi = min(hi(1),domhi(1))
      jhi = min(hi(2),domhi(2))
      khi = min(hi(3),domhi(3))
      
      call filcc (rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif    
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: p        <=  pressure array
c ::: DIMS(p)   => index extent of p array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
c ::: -----------------------------------------------------------

      subroutine FORT_PRESFILL (p,DIMS(p),domlo,domhi,dx,xlo,time,bc)

      implicit none

      integer    DIMDEC(p)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     p(DIMV(p))
      integer    bc(SDIM,2)

      integer    i, j, k
      integer    ilo, ihi, jlo, jhi, klo, khi
      logical    fix_xlo, fix_xhi, fix_ylo, fix_yhi, fix_zlo, fix_zhi
      logical    per_xlo, per_xhi, per_ylo, per_yhi, per_zlo, per_zhi

      fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
      per_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .eq. INT_DIR)
      fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
      per_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .eq. INT_DIR)
      fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
      per_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .eq. INT_DIR)
      fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)
      per_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .eq. INT_DIR)
      fix_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .ne. INT_DIR)
      per_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .eq. INT_DIR)
      fix_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .ne. INT_DIR)
      per_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .eq. INT_DIR)

      ilo = max(ARG_L1(p),domlo(1))
      jlo = max(ARG_L2(p),domlo(2))
      klo = max(ARG_L3(p),domlo(3))
      ihi = min(ARG_H1(p),domhi(1))
      jhi = min(ARG_H2(p),domhi(2))
      khi = min(ARG_H3(p),domhi(3))

c***************
c  SETTING XLO
c***************

      if (fix_xlo) then
         do i = ARG_L1(p), domlo(1)-1
            do k = klo, khi
               do j = jlo,jhi
                  p(i,j,k) = p(ilo,j,k)
               end do 
            end do
       end do

       if (fix_ylo) then
          do i = ARG_L1(p), domlo(1)-1
               do j = ARG_L2(p), domlo(2)-1
                  do k = klo, khi
                     p(i,j,k) = p(ilo,jlo,k)
                  end do
               end do
          end do

          if (fix_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,jlo,klo)
                     end do
                  end do
               end do
          else if (per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,jlo,k)
                     end do
                  end do
               end do
          end if
          if (fix_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,jlo,khi)
                     end do
                  end do
               end do
          else if (per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,jlo,k)
                     end do
                  end do
               end do
          end if
       end if

       if (fix_yhi) then
          do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  do k = klo, khi
                     p(i,j,k) = p(ilo,jhi,k)
                  end do
               end do
          end do
          if (fix_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,jhi,klo)
                     end do
                  end do
               end do
          else if (per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,jhi,k)
                     end do
                  end do
               end do
          end if
          if (fix_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,jhi,khi)
                     end do
                  end do
               end do
          else if (per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,jhi,k)
                     end do
                  end do
               end do
          end if
       end if

       if (fix_zlo) then
          do i = ARG_L1(p), domlo(1)-1
               do j = jlo, jhi
                  do k = ARG_L3(p), domlo(3)-1
                     p(i,j,k) = p(ilo,j,klo)
                  end do
               end do
          end do
            if (per_ylo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,j,klo)
                     end do
                  end do
               end do
            end if
            if (per_yhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,j,klo)
                     end do
                  end do
               end do
            end if

       end if

       if (fix_zhi) then
          do i = ARG_L1(p), domlo(1)-1
               do j = jlo, jhi
                  do k = domhi(3)+1, ARG_H3(p)
                     p(i,j,k) = p(ilo,j,khi)
                  end do
               end do
          end do
            if (per_ylo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,j,khi)
                     end do
                  end do
               end do
            end if
            if (per_yhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,j,khi)
                     end do
                  end do
               end do
            end if
       end if
 
         if (per_ylo) then
               do i = ARG_L1(p), domlo(1)-1
                  do k = klo,khi
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do k = klo,khi
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
 
         if (per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = jlo,jhi
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = jlo,jhi
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_ylo .and. per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
       end if

         if (per_ylo .and. per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
       end if

         if (per_yhi .and. per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
       end if

         if (per_yhi .and. per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
       end if

      end if            

c*****************************************************************************
c SETTING XHI
c*****************************************************************************

      if (fix_xhi) then
         do i = domhi(1)+1, ARG_H1(p)
            do k = klo, khi
               do j = jlo,jhi
                  p(i,j,k) = p(ihi,j,k)
               end do
            end do
       end do

       if (fix_ylo) then
          do i = domhi(1)+1, ARG_H1(p)
               do j = ARG_L2(p), domlo(2)-1
                  do k = klo, khi
                     p(i,j,k) = p(ihi,jlo,k)
                  end do
               end do
          end do

          if (fix_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,jlo,klo)
                     end do
                  end do
               end do
          else if (per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,jlo,k)
                     end do
                  end do
               end do
          end if
          if (fix_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,jlo,khi)
                     end do
                  end do
               end do
          else if (per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,jlo,k)
                     end do
                  end do
               end do
          end if
       end if
       if (fix_yhi) then
          do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  do k = klo, khi
                     p(i,j,k) = p(ihi,jhi,k)
                  end do
               end do
          end do
          if (fix_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,jhi,klo)
                     end do
                  end do
               end do
          else if (per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,jhi,k)
                     end do
                  end do
               end do
          end if
          if (fix_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,jhi,khi)
                     end do
                  end do
               end do
          else if (per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,jhi,k)
                     end do
                  end do
               end do
          end if
       end if

       if (fix_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
               do j = jlo, jhi
                  do k = ARG_L3(p), domlo(3)-1
                     p(i,j,k) = p(ihi,j,klo)
                  end do
               end do
          end do
            if (per_ylo) then
             do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,j,klo)
                     end do
                  end do
               end do
            end if
            if (per_yhi) then
             do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,j,klo)
                     end do
                  end do
               end do
            end if

       end if

       if (fix_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
               do j = jlo, jhi
                  do k = domhi(3)+1, ARG_H3(p)
                     p(i,j,k) = p(ihi,j,khi)
                  end do
               end do
          end do
            if (per_ylo) then
             do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,j,khi)
                     end do
                  end do
               end do
            end if
            if (per_yhi) then
             do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,j,khi)
                     end do
                  end do
               end do
            end if
       end if

         if (per_ylo) then
             do i = domhi(1)+1, ARG_H1(p)
                  do k = klo,khi
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
             do i = domhi(1)+1, ARG_H1(p)
                  do k = klo,khi
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_zlo) then
             do i = domhi(1)+1, ARG_H1(p)
                  do j = jlo,jhi
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
              do i = domhi(1)+1, ARG_H1(p)
                  do j = jlo,jhi
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if


         if (per_ylo .and. per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_ylo .and. per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

      end if            

c*****************************************************************************
c SETTING YLO
c*****************************************************************************

      if (fix_ylo) then
         do j = ARG_L2(p), domlo(2)-1
            do k = klo, khi
               do i = ilo, ihi
                  p(i,j,k) = p(i,jlo,k)
               end do
            end do
       end do

       if (fix_zlo) then
          do j = ARG_L2(p), domlo(2)-1
               do k = ARG_L3(p), domlo(3)-1
                  do i = ilo, ihi
                     p(i,j,k) = p(i,jlo,klo)
                  end do
               end do
          end do
            if (per_xlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jlo,klo)
                     end do
                  end do
               end do
            end if
            if (per_xhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jlo,klo)
                     end do
                  end do
               end do
            end if
       end if

       if (fix_zhi) then
          do j = ARG_L2(p), domlo(2)-1
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ilo, ihi
                     p(i,j,k) = p(i,jlo,khi)
                  end do
               end do
          end do
            if (per_xlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jlo,khi)
                     end do
                  end do
               end do
            end if
            if (per_xhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jlo,khi)
                     end do
                  end do
               end do
            end if
       end if

         if (per_xlo) then
               do j = ARG_L2(p), domlo(2)-1
                  do k = klo,khi
                     do i = ARG_L1(p), domlo(1)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do j = ARG_L2(p), domlo(2)-1
                  do k = klo,khi
                     do i = domhi(1)+1, ARG_H1(p)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_zlo) then
               do j = ARG_L2(p), domlo(2)-1
                  do i = ilo,ihi
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do j = ARG_L2(p), domlo(2)-1
                  do i = ilo,ihi
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if


         if (per_xlo .and. per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

      end if            
 
c*****************************************************************************
c SETTING YHI
c*****************************************************************************

      if (fix_yhi) then
         do j = domhi(2)+1, ARG_H2(p)
            do k = klo, khi
               do i = ilo, ihi
                  p(i,j,k) = p(i,jhi,k)
               end do
            end do
       end do

       if (fix_zlo) then
          do j = domhi(2)+1, ARG_H2(p)
               do k = ARG_L3(p), domlo(3)-1
                  do i = ilo, ihi
                     p(i,j,k) = p(i,jhi,klo)
                  end do
               end do
          end do
            if (per_xlo) then
               do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jhi,klo)
                     end do
                  end do
               end do
            end if
            if (per_xhi) then
               do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jhi,klo)
                     end do
                  end do
               end do
            end if
       end if

       if (fix_zhi) then
          do j = domhi(2)+1, ARG_H2(p)
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ilo, ihi
                     p(i,j,k) = p(i,jhi,khi)
                  end do
               end do
          end do
            if (per_xlo) then
               do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jhi,khi)
                     end do
                  end do
               end do
            end if
            if (per_xhi) then
               do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jhi,khi)
                     end do
                  end do
               end do
            end if
       end if

         if (per_xlo) then
               do j = domhi(2)+1, ARG_H2(p)
                  do k = klo,khi
                     do i = ARG_L1(p), domlo(1)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do j = domhi(2)+1, ARG_H2(p)
                  do k = klo,khi
                     do i = domhi(1)+1, ARG_H1(p)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_zlo) then
               do j = domhi(2)+1, ARG_H2(p)
                  do i = ilo,ihi
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do j = domhi(2)+1, ARG_H2(p)
                  do i = ilo,ihi
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

      end if            

c*****************************************************************************
c SETTING ZLO
c*****************************************************************************

      if (fix_zlo) then
         do k = ARG_L3(p), domlo(3)-1
            do j = jlo, jhi
               do i = ilo, ihi
                  p(i,j,k) = p(i,j,klo)
               end do
            end do
       end do

         if (per_xlo) then
               do k = ARG_L3(p), domlo(3)-1
                  do j = jlo,jhi
                     do i = ARG_L1(p), domlo(1)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do k = ARG_L3(p), domlo(3)-1
                  do j = jlo,jhi
                     do i = domhi(1)+1, ARG_H1(p)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_ylo) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = ilo,ihi
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = ilo,ihi
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_ylo) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = ARG_L1(p), domlo(1)-1
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_yhi) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = ARG_L1(p), domlo(1)-1
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_ylo) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = domhi(1)+1, ARG_H1(p)
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_yhi) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = domhi(1)+1, ARG_H1(p)
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

      end if            

c*****************************************************************************
c SETTING ZHI
c*****************************************************************************

      if (fix_zhi) then
         do k = domhi(3)+1, ARG_H3(p)
            do j = jlo, jhi
               do i = ilo, ihi
                  p(i,j,k) = p(i,j,khi)
               end do
            end do
       end do

         if (per_xlo) then
               do k = domhi(3)+1, ARG_H3(p)
                  do j = jlo,jhi
                     do i = ARG_L1(p), domlo(1)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do k = domhi(3)+1, ARG_H3(p)
                  do j = jlo,jhi
                     do i = domhi(1)+1, ARG_H1(p)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_ylo) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ilo,ihi
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ilo,ihi
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if


         if (per_xlo .and. per_ylo) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ARG_L1(p), domlo(1)-1
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_yhi) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ARG_L1(p), domlo(1)-1
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_ylo) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = domhi(1)+1, ARG_H1(p)
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_yhi) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = domhi(1)+1, ARG_H1(p)
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

      end if            

c*****************************************************************************

      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine FORT_CALCDIVU(divu,DIMS(divu),plo,phi,xlo,dx,time)

c ... This routine is for user defined divu, for debugging purposes
c     only.
      implicit none

      integer    DIMDEC(divu)
      integer    dlo(SDIM), dhi(SDIM)
      integer    plo(SDIM), phi(SDIM)
      REAL_T     xlo(SDIM), dx(SDIM)
      REAL_T     divu(DIMV(divu))

      integer    i,j,k
      REAL_T     x,time

      dlo(1) = divu_l1
      dlo(2) = divu_l2
      dlo(3) = divu_l3
      dhi(1) = divu_h1
      dhi(2) = divu_h2
      dhi(3) = divu_h3

#if 0
      do k=dlo(3),dhi(3)
         do j=dlo(2),dhi(2)
            do i = dlo(1), dhi(1)
               x = xlo(1) + (i-plo(1)+.5)*dx(1)
               if (x.le..005) then
                  divu(i,j,k) = zero
               else
c     divu(i,j,k) = one
c     divu(i,j,k) = min(one,time/.05)
c     divu(i,j,k) = min(two,one+time/.05)
                  divu(i,j,k) = time
               endif
            enddo
         enddo
      enddo
#endif
#if 0
      do k=dlo(3),dhi(3)
         do j=dlo(2),dhi(2)
            do i = dlo(1), dhi(1)
               divu(i,j,k) = time/.0125
            enddo
         enddo
      enddo
#endif
      do k=dlo(3),dhi(3)
         do j=dlo(2),dhi(2)
            do i = dlo(1), dhi(1)
               x = xlo(1) + (i-plo(1)+.5)*dx(1)
               divu(i,j,k) = cvmgt(one,zero,x.GE.0.5)
            enddo
         enddo
      enddo

      end

      subroutine FORT_CHECK_TURB(add_turb_var,turb_length,
     &                           box_offset,turb_scale_var,delta)
      implicit none
#include "probdata.H"
      integer add_turb_var
      REAL_T turb_length(SDIM)
      integer box_offset(SDIM)
      REAL_T delta(SDIM)
      REAL_T turb_scale_var

      add_turb_var = 0

      if (add_turb .and. zturb_top .gt. zero) then

         add_turb_var = 1
         turb_scale_var = turb_scale

         turb_length(1) = domnhi(1)-domnlo(1)
         turb_length(2) = domnhi(2)-domnlo(2)
         turb_length(3) = zturb_top - zturb_bot
         
         box_offset(1) = 0
         box_offset(2) = 0
         box_offset(3) = int(zturb_bot/delta(3))
      endif

      end

      subroutine FORT_RADLOSS(lo,hi,rad,DIMS(rad),
     &                        T,DIMS(T),Y,DIMS(Y),dx,Patm,time)
      implicit none
#include "ckdwrk.H"
#include "probdata.H"
      integer DIMDEC(rad)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer lo(SDIM), hi(SDIM)
      REAL_T  rad(DIMV(rad))
      REAL_T  T(DIMV(T))
      REAL_T  Y(DIMV(Y),1)
      REAL_T  dx(SDIM), Patm, time

      integer i, j, k
      REAL_T xx, yy, zz, r, dr, dz
      
      do k = lo(3),hi(3)
         zz = float(k)*dx(3)+domnlo(3)
         do j = lo(2),hi(2)
            yy = float(j)*dx(2)+domnlo(2)
            do i = lo(1),hi(1)
               xx = float(i)*dx(1)+domnlo(1)

c  Compute the distance to a disk at z=hotspot_z, radius hotspot_r and thickness
c     of hotspot_dz
               dr = SQRT(xx*xx+yy*yy) - hotspot_r
               if (dr .LT. zero) then
                  if (zz.GT.hotspot_z+hotspot_dz) then
                     r = zz - hotspot_z+hotspot_dz
                  else if (zz.LT.hotspot_z) then
                     r = hotspot_z - zz
                  else
                     r = zero
                  end if
               else if ((zz.GT.hotspot_z).AND.(zz.LT.hotspot_z+hotspot_dz)) then
                  r = MAX(dr,zero)
               else
                  dz = MIN(ABS(hotspot_z+hotspot_dz-zz),ABS(hotspot_z-zz))
                  r = SQRT(dr*dr+dz*dz)                  
               end if
               rad(i,j,k) = half*(one-tanh(two*(r-1.5*hotspot_w)/hotspot_w))
     &                                       *hotspot_htc*(hotspot_T - T(i,j,k))
            end do
         end do
      end do
      end

      subroutine FORT_DERVRADIAL (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,
     $                         level,grid_no)
      implicit none
#include "probdata.H"

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k,n
      integer    nxlo, nxhi, nylo, nyhi, nzlo,nzhi
      REAL_T     x,y,r
      
c     ::::: lets punt if not in domain interior
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERVRADIAL: outside domain")
      endif

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            y = (float(j)+0.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x =( float(i)+0.5d0)*delta(1)+domnlo(1)
               r = sqrt(x**2 + y**2)
               e(i,j,k,1) = (x*dat(i,j,k,1) + y*dat(i,j,k,2))/r
            enddo
         enddo
      enddo
      
      end

      subroutine FORT_DERVTHETA (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,
     $                         level,grid_no)
      implicit none
#include "probdata.H"

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k,n
      integer    nxlo, nxhi, nylo, nyhi, nzlo,nzhi
      REAL_T     x,y,r
      
c     ::::: lets punt if not in domain interior
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERVTHETA: outside domain")
      endif

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            y = (float(j)+0.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+0.5d0)*delta(1)+domnlo(1)
               r = sqrt(x**2 + y**2)
               e(i,j,k,1) = (x*dat(i,j,k,2) - y*dat(i,j,k,1))/r
            enddo
         enddo
      enddo

      end

