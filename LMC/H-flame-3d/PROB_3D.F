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
#include "ChemDriver_F.H"
#include "infl_frc.H"

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

#define SDIM 3

#define USE_PMF_CH4
#undef  USE_PMF_CH4

      subroutine FORT_ACTIVECONTROL(coft,time,dt,myproc,step,restart)
      implicit none
      REAL_T coft,time,dt
      integer myproc,step,restart
      end

      subroutine FORT_HACK(lo, hi, dat, DIMS(dat),
     &                     datc, DIMS(datc), nc, mapO2N, lenmap, ratio)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(dat)
      integer DIMDEC(datc)
      integer nc, lenmap, mapO2N(lenmap), ratio
      REAL_T dat(DIMV(dat),nc)
      REAL_T datc(DIMV(datc),nc)
      call bl_abort('Should not in FORT_HACK for this problem...')
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
#include "cdwrk.H"
#include "bc.H"
#include "INFL_FORCE_F.H"

      integer lo(SDIM), hi(SDIM), i, j
      INTEGER dimFile(3)
      REAL_T dxFile(3)
      REAL_T rand

      namelist /fortin/ vorterr, temperr, adverr, tempgrad, 
     &                  flametracval, probtype,
     &			max_temp_lev, max_vort_lev, 
     &                  splitx, splity, splitz, traceSpecVal,
     &                  xfrontw, yfrontw, zfrontw, 
     &                  refine_nozzle, refine_nozzle_x, 
     &                  refine_nozzle_y, refine_nozzle_z, nozzle_width,
     &                  blobx, bloby, blobr, blobz, xcen,
     &                  v_strength,v_width,v_xcen,v_ycen,v_zcen,v_cl_x,
     &                  zturb_top, zturb_bot, turb_scale,add_turb,
     &                  vel_frame, pmf_offset, direction,
     &                  forceInflow, numInflPlanesStore, strmwse_dir,
     $                  forceLo, forceHi, flct_file, 
     $                  T_in, XH2_in, XO2_in, T_out, XH2_out, XO2_out,
     $                  blobrc, blobs, turbr, turbw
      namelist /heattransin/ pamb, dpdt_factor
c
c      Build `probin' filename -- the name of file containing fortin namelist.
c
      integer maxlen, isioproc
      parameter (maxlen=256)

      character probin*(maxlen)
      REAL_T y, yEval, pmf_H_vals(maxspec+3)
      integer n

      call bl_pd_is_ioproc(isioproc)

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
      
      vorterr = 1.e20
      temperr = zero
      adverr = 1.e20
      tempgrad  = 50.0d0
      flametracval = 0.0001d0
      probtype = 1
      max_temp_lev = 0
      max_vort_lev = 0
      splitx = zero
      splity = zero
      splitz = zero
      traceSpecVal = 1.d-14
      xfrontw = fourth*half*(problo(1)+probhi(1))
      yfrontw = zero
      zfrontw = xfrontw
      refine_nozzle = 0
      refine_nozzle_z = zero
      refine_nozzle_y = zero
      refine_nozzle_x = zero
      nozzle_width = -one
      blobx = -one
      bloby = -one
      blobr = -one
      blobz = -one
      xcen = problo(1)
      v_strength = zero
      v_width = zero
      v_xcen = zero
      v_ycen = zero
      v_zcen = zero
      v_cl_x = probhi(1)
      pamb = 101325.d0
      dpdt_factor = 0.3d0
      zturb_top = zero
      zturb_bot = zero
      turb_scale = one
      vel_frame = zero
      pmf_offset = zero
      add_turb = .false.
      direction = "z"
      turbr = one
      turbw = one

      forceInflow = .FALSE.
      numInflPlanesStore = -1
      forceLo = .TRUE.
      forceHi = .FALSE.
      strmwse_dir = FLCT_ZVEL
      flct_file = " "

      T_in = 1000.
      XH2_in  = zero
      XO2_in  = zero
      T_out = 1000.
      XH2_out = zero
      XO2_out = zero
      blobrc  = one
      blobs   = one
      
      read(untin,fortin)
      
      if (max_vort_lev.lt.0) max_vort_lev=max_temp_lev
      
      read(untin,heattransin)

   10 close(unit=untin)

      if ((forceInflow .eqv. .FALSE.).and. .not.((probtype.eq.30).or.(probtype.eq.32))) then
         forceLo = .FALSE.
         forceHi = .FALSE.
      else
         if (flct_file.ne." ") then
#define FF_UNIT 20
            write(6,*) '...initializing turbulence, reading header info'
            open(FF_UNIT, file=trim(flct_file)//'/HDR',form='formatted',status='old')
            call RD_FLCTHD(FF_UNIT, dimFile, probSizeFile, dxFile)
            close(FF_UNIT)
         endif
      endif
c
c     set convection velocity
c
      y = 0.0d0
      yEval = (y-pmf_offset)*1.d2

#ifdef USE_PMF_CH4
      call pmf_CH4(yEval,yEval,pmf_H_vals,n)
#else
      call pmf_H(yEval,pmf_H_vals,n)
#endif

      if (n.ne.Nspec+3) call bl_abort('PROBINIT: n .ne. Nspec+3')

      convVel = pmf_H_vals(2)/1.d2 - vel_frame

c     Load domain dimensions into common, and set up boundary functions
      domnlo(1) = problo(1)
      domnlo(2) = problo(2)
      domnlo(3) = problo(3)
      domnhi(1) = probhi(1)
      domnhi(2) = probhi(2)
      domnhi(3) = probhi(3)
      
      call setupbc_H()
      bcinit = .true.
c      
c     Set random numbers for Ylm perts
c
      do j = lmodemin,lmodemax
         do i = 0,j
            call blutilrand(alphalm(i,j))
            call blutilrand(betalm(i,j))
            call blutilrand(gammalm(i,j))
         end do
      end do

      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
      end if
      
      end

c ::: -----------------------------------------------------------

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
#define ZLO 4
#define ZHI 5

      subroutine setupbc_H()
      implicit none
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "htdata.H"
      
      REAL_T Patm, x, pmf_H_vals(maxspec+3), sum, Y_pmf(maxspec)
      integer i, zone, n, lo(SDIM), hi(SDIM), iN2, getZone
      character*(maxspnml) name
      data lo / 1, 1, 1 /
      data hi / 1, 1, 1 /
      
c     Should we read this data in somehow??
      REAL_T Y_H2(Nzones),Y_O2(Nzones),Y_X(Nzones),Y_N2(Nzones),T(Nzones)
      REAL_T mdot(Nzones),Y_H2O(Nzones), zEval
c     ZONES:               1       2         3       4      5
      data Y_H2  /    0.05d0,   zero,   0.05d0,   zero,   zero /
      data Y_O2  /    zero,    fourth,   zero,   fourth,  zero /
      data Y_H2O /    zero,     zero,    zero,    zero,   half /
      data T     /  1050.d0,   1000.d0, 298.d0, 298.d0,  2.5d3 /
      data mdot  /    two,      four,    one,     one,   .05d0 /
      
      Patm = pamb / 101325.0d0

      if (probtype .eq. 8 .or. probtype .eq. 28 .or. probtype .eq. 29) then
c     For probtype == 8, set boundary conditions from premixed flame profile 
c     z direction
         x = domnlo(1)
         zone = getZone(x,domnlo(2),domnlo(3))
         iN2 = -1
         do n = 1,Nspec
            call get_spec_name(name,n)
            if (name .eq. 'N2' ) iN2 = n
         end do

         zEval = (x-pmf_offset)*1.d2

#ifdef USE_PMF_CH4
         call pmf_CH4(zEval,zEval,pmf_H_vals,n)
#else
         call pmf_H(zEval,pmf_H_vals,n)
#endif

         if (n.ne.Nspec+3) call bl_abort('SETUPBC: n .ne. Nspec+3')

         call FORT_MOLETOMASS(lo,hi,
     &        pmf_H_vals(4), ARLIM(lo), ARLIM(hi),
     &        Y_pmf,         ARLIM(lo), ARLIM(hi))
         
         if (probtype .eq. 8 .or. probtype .eq. 28) then
            do n=0,Nspec-1
               Y_bc(n,zone) = Y_pmf(n+1)
            end do
         else if (probtype .eq. 29) then
            do n=0,Nspec-1
               Y_bc(n,zone) = zero
            enddo
            Y_bc(iN2,zone) = one
         endif
         if (iN2.gt.0) then
            sum = zero
            do n = 1,Nspec
               if (n.ne.iN2) sum = sum+Y_bc(n-1,zone)
            end do
            Y_bc(iN2-1,zone) = one - sum
         end if
         T_bc(zone) = pmf_H_vals(1)
         u_bc(zone) = zero
         v_bc(zone) = zero
         w_bc(zone) = pmf_H_vals(2) / 1.d2 - vel_frame
         call FORT_RHOfromPTY(lo,hi,
     &                        rho_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        T_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        Y_bc(0,zone),ARLIM(lo),ARLIM(hi), Patm)
         call FORT_HMIXfromTY(lo,hi,
     &                        h_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        T_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        Y_bc(0,zone),ARLIM(lo),ARLIM(hi))
         
         
      else

         do zone = 1,Nzones
            Y_X(zone) = traceSpecVal
            if (Y_H2(zone).eq.zero) Y_H2(zone)=traceSpecVal
            if (Y_O2(zone).eq.zero) Y_O2(zone)=traceSpecVal
            if (Y_H2O(zone).eq.zero) Y_H2O(zone)=traceSpecVal
            Y_N2(zone) = one - Y_H2(zone) - Y_O2(zone) - Y_H2O(zone)
     &           - (Nspec-4)*Y_X(zone)
         end do
         
c     Initialize zone values
         do zone = 1,Nzones
            do n = 0,Nspec-1
               Y_bc(n,zone) = Y_X(zone)
               call get_spec_name(name, n+1)
               if (name .EQ. 'H2')  Y_bc(n,zone) = Y_H2(zone)
               if (name .EQ. 'O2')  Y_bc(n,zone) = Y_O2(zone)
               if (name .EQ. 'H2O') Y_bc(n,zone) = Y_H2O(zone)
               if (name .EQ. 'N2')  Y_bc(n,zone) = Y_N2(zone)
            end do
            T_bc(zone) = T(zone)
            
            call FORT_RHOfromPTY(lo,hi,
     &           rho_bc(zone),ARLIM(lo),ARLIM(hi),
     &           T_bc(zone),ARLIM(lo),ARLIM(hi),
     &           Y_bc(0,zone),ARLIM(lo),ARLIM(hi), Patm)
            call FORT_HMIXfromTY(lo,hi,
     &           h_bc(zone),ARLIM(lo),ARLIM(hi),
     &           T_bc(zone),ARLIM(lo),ARLIM(hi),
     &           Y_bc(0,zone),ARLIM(lo),ARLIM(hi))
            if (probtype .ne. 1) then
               u_bc(zone) = zero
               v_bc(zone) = zero
               w_bc(zone) = mdot(zone)
            elseif (direction .eq. "z") then
               u_bc(zone) = zero
               v_bc(zone) = zero
               w_bc(zone) = mdot(zone)
            elseif (direction .eq. "y") then
               u_bc(zone) = zero
               v_bc(zone) = mdot(zone)
               w_bc(zone) = zero
            elseif (direction .eq. "x") then
               u_bc(zone) = mdot(zone)
               v_bc(zone) = zero
               w_bc(zone) = zero
            endif
         end do
      end if

      end
         
c ::: -----------------------------------------------------------

      integer function getZone(x, y, z)
      implicit none
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      REAL_T x, y, z
c     
c     probtype:
c     4  ->  Nozzle centered on splitx (zone 2 inside nozzle, else zone 1)
c     5  ->  Nozzle centered at origin ( -same- )
c     else 
c       zone 1: x<splitx, zone 2: x>=splitx
c           
      getZone = 1
      if (probtype .eq. 4) then
         if (ABS(x-splitx) .ge. half*nozzle_width) getZone = 2
      else if (probtype .eq. 5) then
         if (x .ge. half*nozzle_width) getZone = 2
      else if (probtype .eq. 6) then
         getZone = 3
         if (ABS(x-xcen) .ge. splitx) getZone = 4
      else if( (probtype .eq. 2) .or.
     &         (probtype .eq. 3) .or.
     &         (probtype .eq. 7) ) then
         if (x .ge. splitx) getZone = 2
      else if (probtype .eq. 8 .or. probtype .eq. 28
     &        .or. probtype .eq. 29) then
         getZone = 1
      else if (probtype .eq. 1) then
         if (direction .eq. "z" .and. z .ge. splitx) getZone = 2
         if (direction .eq. "y" .and. y .ge. splitx) getZone = 2
         if (direction .eq. "x" .and. x .ge. splitx) getZone = 2
      end if
      end

      subroutine bcfunction(orient, x, y, z, time, u, v, w, rho, Yl, T, h)

      implicit none

      integer orient
      REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h

#include "htdata.H"
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"

      REAL_T Patm, eta, zblend, mdotL, mdotR, mdotb
      integer n, getZone, zone, lo(SDIM), hi(SDIM), zL, zR
      data lo / 1, 1, 1 /
      data hi / 1, 1, 1 /
c
c     A handy statement function: here blend goes from 0 to 1 in x at
c     splitx, over a width of xfrontw.
c
      zblend(x) = half*(one+TANH((x-splitx)/(fourth*xfrontw)))
c     
c     Some boundary conditions may not be not time-dependent, but rather
c     just constant values for each orientation.  We can precompute
c     then (in setupbc), here just retrieve the appropriate values.
c     However, we do not know the size of the domain here, so we cannot
c     call setupbc values ourselves.  Need to assume someone else did it.
c
      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if
      
      if ((probtype .eq. 3).or.(probtype .eq. 7)) then
         zL = 1
         zR = 2
         if (probtype.eq.7) then
            zL = 3
            zR = 4
         end if 
         eta = zblend(ABS(x-xcen))
         T = (one-eta)*T_bc(zL) + eta*T_bc(zR)
         do n = 0,Nspec-1
            Yl(n) = (one-eta)*Y_bc(n,zL) +eta*Y_bc(n,zR)
         enddo
         Patm = pamb / 101325.0d0
         call FORT_RHOfromPTY(lo,hi,rho,ARLIM(lo),ARLIM(hi),
     &                        T,ARLIM(lo),ARLIM(hi),
     &                        Yl,ARLIM(lo),ARLIM(hi), Patm)
         call FORT_HMIXfromTY(lo,hi,h,ARLIM(lo),ARLIM(hi),
     &                        T,ARLIM(lo),ARLIM(hi),
     &                        Yl,ARLIM(lo),ARLIM(hi))

         mdotL = v_bc(zL)*rho_bc(zL)
         mdotR = v_bc(zR)*rho_bc(zR)
         mdotb = (one-eta)*mdotL + eta*mdotR
         u = zero
         v = zero
         w = mdotb / rho
         if (probtype.eq.7) then
            u = zero
            v = zero
            w = (one-eta)*w_bc(zL) + eta*w_bc(zR)
         end if
      else
         if (orient .le. ZHI) then
            zone = getZone(x,y,z)
            u = u_bc(zone)
            v = v_bc(zone)
            w = w_bc(zone)
            rho = rho_bc(zone)
            do n = 0,Nspec-1
               Yl(n) = Y_bc(n,zone)
            end do
            T = T_bc(zone)
            h = h_bc(zone)
         else
            write(6,*) 'No boundary condition for orientation = ', orient
            call bl_abort(" ")
         end if
      end if
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
     &	 	               vel,scal,DIMS(state),press,DIMS(press),
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
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "INFL_FORCE_F.H"

      integer i, j, k, n
      REAL_T  x, y, z, r
      integer slo(SDIM), shi(SDIM)
      character*(maxspnml) name
      REAL_T Patm
      integer spec, getZone, zone
      
      integer zL, zR, zH
      REAL_T eta, xi, xblend, mdotb, mdotL, mdotR, zblend
c      REAL_T r1sq, r2sq
      REAL_T pmf_H_vals(maxspec+3), sum, Y_pmf(maxspec), zEval, zlo, zhi
      REAL_T Xini(maxspec), Yini(maxspec)
      integer iN2, iH2, iO2
      integer statelo1,statelo2,statelo3
      integer lo_chem(SDIM), hi_chem(SDIM)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/
      
      integer DIMDEC(uflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T uflct(:,:,:)
      allocatable uflct
      REAL_T fileDataZLO, fileDataZHI, gridLo(SDIM), convVelSAVE
      integer bc(SDIM,2)
      
      REAL_T theta, phi, rfront, Ylm
      integer l, m, ctr

c     A handy statement function: here blend goes from 0 to 1 in x at
c     splitx, over a width of xfrontw.  The interface is a jump below
c     splitz, and smoothed above it, transitioning over a width zfrontw
      xblend(x,y,z) = half*(one+TANH((x-splitx)/(eighth*xfrontw*
     &                     (one+TANH((z-splitz)/(fourth*zfrontw))))))
c     A handy statement function: here blend goes from 0 to 1 in x at
c     splitx, over a width of xfrontw.  
      zblend(x) = half*(one+TANH((x-splitx)/(fourth*xfrontw)))
c     
c     IC/BCs:
c     1. Constant value, zone = 1 - use direction to specify direction
c     2. Discontinuous at the inflow boundary, transitioning to a
c           widened hyperbolic tangent (zones = 1,2)
c     3. Hyperbolic tangent in x (zones = 1,2), BCs are discontinuous
c     4. Discontinuous across zones (zone=1 inside nozzle centered
c           at split x, width = nozzle, zone=2 otherwise)
c     5. Like 4, but intended for r-z with nozzle centered at origin
c     6. ICs with hyperbolic tangent split between zones 3 and 4 at
c           the inflow, but with a zone 5 blob superimposed with
c           smoothed edges.  BCs are discontinuous.
c     7. Like 6, but smoothed BCs as well
c
c     8. This is a 1D flame (in the z-direction) with a vortex
c
c     28.  This is a 1D flame (in the z-direction) with a patch of
c         turbulence at the bottom of the domain
c
c     29. This is probtype 8, with only N2 present (mass fraction of 1).
c         Everything else has a mass fraction of 0.
c     
c     30. This is a Hydrogen flame ball with outflow everytwhere
c     31. Same as 30, but perturb the surface of the ball, and put stationary field
c     32. Same as 30, but use file giving post-pmf burn profiles
c      
      if ((Temp.gt.0).neqv.(RhoH.gt.0)) then
         call bl_abort('Need both Temp and RhoH, or neither')
      end if
      
      if ((Temp .LT. 0) .OR. (RhoH .LT. 0)) then
         call bl_abort('No IC''s for system without T, RhoH')
      endif
      
      statelo1 = ARG_L1(state)
      statelo2 = ARG_L2(state)
      statelo3 = ARG_L3(state)

      if ( (probtype .eq. 1) .or. (probtype .eq. 4) .or.
     &     (probtype .eq. 5)) then
c
c     Set IC's to constants in each zone
c
         do k = lo(3), hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  zone = getZone(x,y,z)
                  scal(i,j,k,Trac) = zero
                  scal(i,j,k,Density) = rho_bc(zone)
                  scal(i,j,k,Temp) = T_bc(zone)
                  scal(i,j,k,RhoH) = rho_bc(zone)*h_bc(zone)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) = rho_bc(zone)*Y_bc(n,zone)
                  enddo
                  vel(i,j,k,1) = u_bc(zone)
                  vel(i,j,k,2) = v_bc(zone)
                  vel(i,j,k,3) = w_bc(zone)
               enddo
            enddo
         enddo
         
      else if (probtype .eq. 2) then
c     
c     Set IC's to join zones 1 and 2 across an x-interface, making it
c     discontinuous at the zlo boundary, and smoothing across a progressively
c     wider zone with z
c
         zL = 1
         zR = 2
         do k = lo(3), hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  eta = xblend(x,y,z)
                  scal(i,j,k,Trac) = zero
                  scal(i,j,k,Temp) = (one-eta)*T_bc(zL) + eta*T_bc(zR)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) = (one-eta)*Y_bc(n,zL) 
     &                                          +eta*Y_bc(n,zR)
                  enddo
               enddo
            enddo
         enddo
         Patm = pamb / 101325.0d0
         call FORT_RHOfromPTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,Density),  DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state),
     &        Patm)
         call FORT_HMIXfromTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,RhoH),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state))
         mdotL = w_bc(zL)*rho_bc(zL)
         mdotR = w_bc(zR)*rho_bc(zR)
         do k = lo(3), hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  eta = xblend(x,y,z)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) = scal(i,j,k,firstSpec+n)
     &                    *scal(i,j,k,Density)
                  enddo
                  scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
                  mdotb = (one-eta)*mdotL + eta*mdotR
                  vel(i,j,k,1) = zero
                  vel(i,j,k,2) = zero
                  vel(i,j,k,3) = mdotb / scal(i,j,k,Density)
               enddo
            enddo
         enddo
      else if (probtype .eq. 3) then         
c
c     Set IC's to join zones 1 and 2 across an x-interface
c      
         zL = 1
         zR = 2
         do k = lo(3), hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  eta = zblend(x)
                  scal(i,j,k,Trac) = zero
                  scal(i,j,k,Temp) = (one-eta)*T_bc(zL) + eta*T_bc(zR)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) = (one-eta)*Y_bc(n,zL) 
     &                                       +eta*Y_bc(n,zR)
                  enddo
               enddo
            enddo
         enddo
         Patm = pamb / 101325.0d0
         call FORT_RHOfromPTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,Density),  DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state),
     &        Patm)
         call FORT_HMIXfromTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,RhoH),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state))
         mdotL = w_bc(zL)*rho_bc(zL)
         mdotR = w_bc(zR)*rho_bc(zR)
         do k = lo(3), hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  eta = zblend(x)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) = scal(i,j,k,firstSpec+n)
     &                                        *scal(i,j,k,Density)
                  enddo
                  scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
                  mdotb = (one-eta)*mdotL + eta*mdotR
                  vel(i,j,k,1) = zero
                  vel(i,j,k,2) = zero
                  vel(i,j,k,3) = mdotb / scal(i,j,k,Density)
               enddo
            enddo
         enddo
      else if (probtype .eq. 6) then         
c
c     A hot spot of products above a column of fuel extending from the
c       nozzle
c     
c     xcen      center of H2 stream
c     splitx    half-width of H2 stream
c     xfrontw   width of transitions
c     blobx     x of center of hot blob
c     bloby     y of center of hot blob
c     blobr     radius of hot blob
c         
c     ZONE 3: H2 + N2    inside nozzle jet, into blob
c     ZONE 4: O2 + N2    outside of nozzle jet, and above blob
c     ZONE 5  H2O + N2   in hot blob
c      
         zR = 4
         zH = 5
         do k = lo(3), hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               zL = 3
               if (z.gt.blobz) zL = zR
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  eta = zblend(ABS(x-xcen))
                  scal(i,j,k,Trac) = zero
                  scal(i,j,k,Temp) = (one-eta)*T_bc(zL) + eta*T_bc(zR)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) = (one-eta)*Y_bc(n,zL)
     &                                         +eta*Y_bc(n,zR)
                  enddo
c     Superimpose a spot of hot products (values in zone zH)
                  xi = half*(one+TANH((SQRT((x-blobx)**2
     &                                  +(z-blobz)**2)-blobr)/(.25*xfrontw)))
                  scal(i,j,k,Temp) = (one-xi)*T_bc(zH) + xi*scal(i,j,k,Temp)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) =
     &                    (one-xi)*Y_bc(n,zH)+xi*scal(i,j,k,firstSpec+n)
                  enddo
               enddo
            enddo
         enddo
         Patm = pamb / 101325.0d0
         call FORT_RHOfromPTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,Density),  DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state),
     &        Patm)
         call FORT_HMIXfromTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,RhoH),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state))
         do k = lo(3), hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               zL = 3
               if (z.gt.blobz) zL = zR
               mdotL = w_bc(zL)*rho_bc(zL)
               mdotR = w_bc(zR)*rho_bc(zR)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  eta = zblend(ABS(x-xcen))
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) = scal(i,j,k,firstSpec+n)
     &                                        *scal(i,j,k,Density)
                  enddo
                  scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
                  mdotb = (one-eta)*mdotL + eta*mdotR
                  xi = half*(one+TANH((SQRT(x**2 +(z-blobz)**2)
     &                                 -2*splitx)/(.25*xfrontw)))
                  mdotb = (one-xi)*mdotb + w_bc(zH)*rho_bc(zH)
                  vel(i,j,k,1) = zero
                  vel(i,j,k,2) = zero
                  vel(i,j,k,3) = mdotb / scal(i,j,k,Density)
               enddo
            enddo
         enddo
      else if (probtype .eq. 7) then         
         zR = 4
         zH = 5
         do k = lo(3), hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               zL = 3
               if (z.gt.blobz) zL = zR
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  eta = zblend(ABS(x-xcen))
                  scal(i,j,k,Trac) = zero
                  scal(i,j,k,Temp) = (one-eta)*T_bc(zL) + eta*T_bc(zR)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) = (one-eta)*Y_bc(n,zL)
     &                                          +eta*Y_bc(n,zR)
                  enddo
c     Superimpose a spot of hot products (values in zone zH)
c                  xi = half*(one+TANH((SQRT((x-blobx)**2
c     &                 +(z-blobz)**2)-blobr)/(.25*xfrontw)))
                  xi = half*(one+TANH((SQRT(
     &                 (z-blobz)**2)-blobr)/(.25*xfrontw)))
                  scal(i,j,k,Temp) = (one-xi)*T_bc(zH) + xi*scal(i,j,k,Temp)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) =
     &                    (one-xi)*Y_bc(n,zH)+xi*scal(i,j,k,firstSpec+n)
                  enddo
                  vel(i,j,k,1) = zero
                  vel(i,j,k,2) = zero
                  vel(i,j,k,3) = (one-eta)*w_bc(zL) + eta*w_bc(zR)
               enddo
            enddo
         enddo
         Patm = pamb / 101325.0d0
         call FORT_RHOfromPTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,Density),  DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state),
     &        Patm)
         call FORT_HMIXfromTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,RhoH),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state))
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) = scal(i,j,k,firstSpec+n)
     &                                        *scal(i,j,k,Density)
                  enddo
                  scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
               enddo
            enddo
         enddo
      else if (probtype .eq. 8) then
         iN2 = -1
         do n = 1,Nspec
            call get_spec_name(name,n)
            if (name .eq. 'N2' ) iN2 = n
         end do
c     no vorticity
c     Set vorticity into pressure field
c     do j = ARG_L2(press),ARG_H2(press)
c            y = float(j)*delta(2)+domnlo(2)
c           do i = ARG_L1(press),ARG_H1(press)
c            x = float(i)*delta(1)+domnlo(1)
c            r1sq = (x-v_xcen)**2            + (y-v_ycen)**2
c               r2sq = (x-(two*v_cl_x-v_xcen))**2 + (y-v_ycen)**2
c               press(i,j) = v_strength*
c     &              (dexp(-v_width*r1sq) - dexp(-v_width*r2sq)) 
c            end do
c     end do
         
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

         loFlctArray(strmwse_dir) = 1
         hiFlctArray(strmwse_dir) = 1
         call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
         allocate(uflct(DIMV(uflct)))

         do i = 1, SDIM
            bc(i,1) = INT_DIR
            bc(i,2) = INT_DIR
         end do

         do k= lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            if (flct_file.ne." ") then           
               convVelSAVE = convVel
               convVel = one
               call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo,
     &              delta, z - fileDataZLO, bc, domnlo, domnhi)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     vel(i,j,k,1) = uflct(i,j,1)
                  end do
               end do
               call INFL_FILL(FLCT_YVEL, DIMS(uflct), uflct, xlo,
     &              delta, z - fileDataZLO, bc, domnlo, domnhi)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     vel(i,j,k,2) = uflct(i,j,1)
                  end do
               end do
               call INFL_FILL(FLCT_ZVEL, DIMS(uflct), uflct, xlo,
     &              delta, z - fileDataZLO, bc, domnlo, domnhi)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     vel(i,j,k,3) = uflct(i,j,1)
                  end do
               end do
               convVel = convVelSAVE
            else
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     do n=1,SDIM
                        vel(i,j,k,n) = zero
                     end do
                  end do
               end do
            end if

            zEval = (z-pmf_offset)*1.d2

#ifdef USE_PMF_CH4
            zlo = (    k*delta(3) + domnlo(3) - pmf_offset)*1.d2
            zhi = ((k+1)*delta(3) + domnlo(3) - pmf_offset)*1.d2
            call pmf_CH4(zlo,zhi,pmf_H_vals,n)
#else
            call pmf_H(zEval,pmf_H_vals,n)
#endif

            if (n.ne.Nspec+3) call bl_abort('INITDATA: n .ne. Nspec+3')

            call FORT_MOLETOMASS(lo_chem,hi_chem,
     &           pmf_H_vals(4), ARLIM(lo_chem), ARLIM(hi_chem),
     &           Y_pmf,         ARLIM(lo_chem), ARLIM(hi_chem))
            
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  scal(i,j,k,Trac) = zero
c                  scal(i,j,Trac) = fourth*(press(i,j  ) + press(i+1,j)
c     &                 +                   press(i,j+1) + press(i+1,j+1))
                  scal(i,j,k,Temp) = pmf_H_vals(1)
                  do n = 1,Nspec
                     scal(i,j,k,firstSpec+n-1) = Y_pmf(n)
                  end do
                  if (iN2.gt.0) then
                     sum = zero
                     do n = 1,Nspec
                        if (n.ne.iN2) sum = sum+scal(i,j,k,firstSpec+n-1)
                     end do
                     scal(i,j,k,firstSpec+iN2-1) = one - sum
                  end if
                  vel(i,j,k,3) = vel(i,j,k,3) + 
     &                   pmf_H_vals(2) / 1.d2 - vel_frame
               enddo
            enddo
         enddo
         Patm = pamb / 101325.0d0
         call FORT_RHOfromPTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,Density),  DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state),
     &        Patm)
         call FORT_HMIXfromTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,RhoH),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state))
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) 
     &                    = scal(i,j,k,firstSpec+n)*scal(i,j,k,Density)
                  enddo
                  scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
               enddo
            enddo
         enddo

      if (allocated(uflct)) deallocate(uflct)

      else if (probtype .eq. 28 .or. probtype .eq. 29) then
         iN2 = -1
         do n = 1,Nspec
            call get_spec_name(name,n)
            if (name .eq. 'N2' ) iN2 = n
         end do
         
         do k= lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            
            zEval = (z-pmf_offset)*1.d2  

#ifdef USE_PMF_CH4
            zlo = (    k*delta(3) + domnlo(3) - pmf_offset)*1.d2
            zhi = ((k+1)*delta(3) + domnlo(3) - pmf_offset)*1.d2
            call pmf_CH4(zlo,zhi,pmf_H_vals,n)
#else
            call pmf_H(zEval,pmf_H_vals,n)
#endif

            if (n.ne.Nspec+3) call bl_abort('INITDATA: n .ne. Nspec+3')

            call FORT_MOLETOMASS(lo_chem,hi_chem,
     &           pmf_H_vals(4), ARLIM(lo_chem), ARLIM(hi_chem),
     &           Y_pmf,         ARLIM(lo_chem), ARLIM(hi_chem))

            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  scal(i,j,k,Trac) = zero
c                  scal(i,j,Trac) = fourth*(press(i,j  ) + press(i+1,j)
c     &                 +                   press(i,j+1) + press(i+1,j+1))
                  scal(i,j,k,Temp) = pmf_H_vals(1)

                  if (probtype .eq. 28) then
                     do n = 1,Nspec
                        scal(i,j,k,firstSpec+n-1) = Y_pmf(n)
                     end do
                  else if (probtype .eq. 29) then
                     do n = 1,NSpec
                        scal(i,j,k,firstSpec+n-1) = zero
                     enddo
                     scal(i,j,k,firstSpec+iN2-1) = one
                  endif
                  if (iN2.gt.0) then
                     sum = zero
                     do n = 1,Nspec
                        if (n.ne.iN2) sum = sum+scal(i,j,k,firstSpec+n-1)
                     end do
                     scal(i,j,k,firstSpec+iN2-1) = one - sum
                  end if
                  vel(i,j,k,1) = zero
                  vel(i,j,k,2) = zero
                  vel(i,j,k,3) = pmf_H_vals(2) / 1.d2
               enddo
            enddo
         enddo

         Patm = pamb / 101325.0d0
         call FORT_RHOfromPTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,Density),  DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state),
     &        Patm)
         call FORT_HMIXfromTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,RhoH),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state))
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) 
     &                    = scal(i,j,k,firstSpec+n)*scal(i,j,k,Density)
                  enddo
                  scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
               enddo
            enddo
         enddo
      else if ((probtype .eq. 30).or.(probtype .eq. 32)) then
         iN2 = -1
         do n = 1,Nspec
            call get_spec_name(name,n)
            if (name .eq. 'N2' ) iN2 = n
            if (name .eq. 'H2' ) iH2 = n
            if (name .eq. 'O2' ) iO2 = n
         end do

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
         convVel = one
         strmwse_dir = 3
         loFlctArray(strmwse_dir) = 1
         hiFlctArray(strmwse_dir) = 1
         call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
         allocate(uflct(DIMV(uflct)))

c     HACK!!    make up some boundary conditions
         do i = 1, SDIM
            bc(i,1) = EXT_DIR
            bc(i,2) = EXT_DIR
         end do
            
         do k= lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n=1,SDIM
                     vel(i,j,k,n) = zero
                  end do
               end do
            end do
            
            z = (float(k)+.5)*delta(3)+domnlo(3)
            
            if (flct_file.ne." ") then

               if ((z.GE.fileDataZLO).AND.(z.LE.fileDataZHI)) then
                  
                  call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo,
     &                 delta, z - fileDataZLO, bc, domnlo, domnhi)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        vel(i,j,k,1) = ten*uflct(i,j,1)
                     end do
                  end do
                  
                  call INFL_FILL(FLCT_YVEL, DIMS(uflct), uflct, xlo,
     &                 delta, z - fileDataZLO, bc, domnlo, domnhi)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        vel(i,j,k,2) = ten*uflct(i,j,1)
                     end do
                  end do
                  
                  call INFL_FILL(FLCT_ZVEL, DIMS(uflct), uflct, xlo,
     &                 delta, z - fileDataZLO, bc, domnlo, domnhi)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        vel(i,j,k,3) = ten*uflct(i,j,1)
                     end do
                  end do
               end if
            end if
            
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  scal(i,j,k,Trac) = zero
                  
                  r = SQRT((x-blobx)**2+(y-bloby)**2+(z-blobz)**2)
                  
                  if (probtype.eq.30) then
                     
                     eta = half*(one+TANH((r-blobrc)*Pi*blobs/blobr))                 
                     
                     scal(i,j,k,Temp) = eta*T_out + (one-eta)*T_in
                     if (iH2.gt.0) Xini(iH2) = eta*XH2_out + (one-eta)*XH2_in
                     if (iO2.gt.0) Xini(iO2) = eta*XO2_out + (one-eta)*XO2_in
                     
                  else
                     
                     call ball(r*1.d2,r*1.d2,pmf_H_vals,n)
                     scal(i,j,k,Temp) = pmf_H_vals(2)
                     do n = 1,Nspec
                        Xini(n) = pmf_H_vals(3+n)
                     end do
                     
                  end if
                  
                  eta = half*(one-TANH(two*(r-turbr)/turbw))
                  do n = 1,SDIM
                     vel(i,j,k,n) = eta*vel(i,j,k,n)
                  end do
                     
                  if (iN2.gt.0) then
                     sum = zero
                     do n = 1,Nspec
                        if (n.ne.iN2) sum = sum+Xini(n)
                     end do
                     Xini(iN2) = one - sum
                  end if
                  
                  call FORT_MOLETOMASS(lo_chem,hi_chem,
     &                    Xini, ARLIM(lo_chem), ARLIM(hi_chem),
     &                    Yini, ARLIM(lo_chem), ARLIM(hi_chem))
                     
                  do n = 1,Nspec
                     scal(i,j,k,firstSpec+n-1) = Yini(n)
                  end do

               end do
            end do
         end do
         
         if (allocated(uflct)) deallocate(uflct)
         
         Patm = pamb / 101325.0d0
         call FORT_RHOfromPTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,Density),  DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state),
     &        Patm)
         call FORT_HMIXfromTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,RhoH),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state))
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) 
     &                    = scal(i,j,k,firstSpec+n)*scal(i,j,k,Density)
                  enddo
                  scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
               enddo
            enddo
         enddo
      else if (probtype .eq. 31) then
         iN2 = -1
         do n = 1,Nspec
            call get_spec_name(name,n)
            if (name .eq. 'N2' ) iN2 = n
            if (name .eq. 'H2' ) iH2 = n
            if (name .eq. 'O2' ) iO2 = n
         end do

         do k= lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n=1,SDIM
                     vel(i,j,k,n) = zero
                  end do
               end do
            end do
         end do
         
         do k = lo(3), hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  scal(i,j,k,Trac) = zero

                  r = SQRT((x-blobx)**2+(y-bloby)**2+(z-blobz)**2)
                  if (r.LE.delta(1)) then
                     eta = zero
                  else
                     theta = ATAN2(SQRT((x-blobx)**2+(y-bloby)**2)/r,(z-blobz)/r)-Pi
                     phi = ATAN2((y-bloby)/r,(x-blobx)/r)
                     
                     rfront = zero
                     ctr = 0
                     do l = lmodemin,lmodemax
                        do m = 0,l
                           ctr = ctr + 1
                           rfront = rfront + 
     &                          betalm(l,m)*EXP(-float(l-lpeak)**2/(.15*lpeak)**2)
     &                          * Ylm(l,m,phi+alphalm(l,m)*2*Pi,
     &                                COS(theta+gammalm(l,m)*2*Pi))
                        end do
                     end do
                     if (ctr.GT.0) rfront = rfront/ctr
                        
                     rfront = blobrc*(one + rfront)
                     
                     eta = half*(one+TANH((r-rfront)*Pi*blobs/blobr))                 
                  endif
                  
                  
                  scal(i,j,k,Temp) = eta*T_out + (one-eta)*T_in
                  if (iH2.gt.0) Xini(iH2) = eta*XH2_out + (one-eta)*XH2_in
                  if (iO2.gt.0) Xini(iO2) = eta*XO2_out + (one-eta)*XO2_in
                  if (iN2.gt.0) then
                     sum = zero
                     do n = 1,Nspec
                        if (n.ne.iN2) sum = sum+Xini(n)
                     end do
                     Xini(iN2) = one - sum
                  end if
                  
                  call FORT_MOLETOMASS(lo_chem,hi_chem,
     &                 Xini, ARLIM(lo_chem), ARLIM(hi_chem),
     &                 Yini, ARLIM(lo_chem), ARLIM(hi_chem))
                  
                  do n = 1,Nspec
                     scal(i,j,k,firstSpec+n-1) = Yini(n)
                  end do
                  
               end do
            end do
         end do
         
         Patm = pamb / 101325.0d0
         call FORT_RHOfromPTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,Density),  DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state),
     &        Patm)
         call FORT_HMIXfromTY(lo,hi,
     &        scal(statelo1,statelo2,statelo3,RhoH),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,Temp),     DIMS(state),
     &        scal(statelo1,statelo2,statelo3,firstSpec),DIMS(state))
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n = 0,Nspec-1
                     scal(i,j,k,firstSpec+n) 
     &                    = scal(i,j,k,firstSpec+n)*scal(i,j,k,Density)
                  enddo
                  scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
               enddo
            enddo
         enddo
      end if
      
      end
      
      REAL_T function Ylm(l,m,phi,costheta)
      REAL_T phi, costheta, facN, facD, plgndr
      integer m, l, n
      facN = l - m
      if (facN .EQ. 0) facN=1
      do n = l-m-1,2,-1
         facN = facN*n
      end do
      facD = l+m
      if (facD.EQ.0) facD=1
      do n = l+m-1,2,-1
         facD = facD*n
      end do
      Ylm = SQRT((2*l+1)*facN/(4*Pi*facD))*plgndr(l,m,costheta)*COS(m*phi)
      end
      
      REAL_T FUNCTION plgndr(l,m,x)
      INTEGER l,m
      REAL_T x
      INTEGER i,ll
      REAL_T fact,pll,pmm,pmmp1,somx2
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) then
         call bl_abort('bad arguments in plgndr')
      endif
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
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
      REAL_T  diff(DIMV(diff),1)
      REAL_T  dx(SDIM)
      REAL_T  problo(SDIM)
      
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

      call bl_abort('DENERROR: should no be here')
      
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
      REAL_T    x, y, z
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

      if (refine_nozzle .eq. 1) then
         do k = lo(3), hi(3)
            z = (float(k)+.5)*dx(3)+problo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*dx(2)+problo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*dx(1)+problo(1)
                  in_refine_zone = (y - problo(2)) .le. refine_nozzle_y
                  in_refine_zone = in_refine_zone .and.
     &                  (z - problo(3)) .le. refine_nozzle_z
                  in_refine_zone = in_refine_zone .and. 
     &                  (ABS(x-xcen) .le. refine_nozzle_x)
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
c
c     ::::: refine where there is temperature gradient
c
      if (level .lt. max_temp_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = abs(temperature(i+1,j,k,1) - temperature(i-1,j,k,1))
                  ay = abs(temperature(i,j+1,k,1) - temperature(i,j-1,k,1))
                  az = abs(temperature(i,j,k+1,1) - temperature(i,j,k-1,1))
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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(den)
      hi(1) = ARG_H1(den)
      lo(2) = ARG_L2(den)
      hi(2) = ARG_H2(den)
      lo(3) = ARG_L3(den)
      hi(3) = ARG_H3(den)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      klo = max(lo(3),domlo(3))
      khi = min(hi(3),domhi(3))
      
      call filcc (den,DIMS(den),domlo,domhi,delta,xlo,bc)
      
c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(XHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(YLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(YHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(ZLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(ZHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(temp)
      hi(1) = ARG_H1(temp)
      lo(2) = ARG_L2(temp)
      hi(2) = ARG_H2(temp)
      lo(3) = ARG_L3(temp)
      hi(3) = ARG_H3(temp)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      klo = max(lo(3),domlo(3))
      khi = min(hi(3),domhi(3))
      
      call filcc (temp,DIMS(temp),domlo,domhi,delta,xlo,bc)
      
c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(XHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(YLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(YHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(ZLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(ZHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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
     &                         xlo,time,bc,id)
      implicit none
      integer DIMDEC(rhoh), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM), id
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoh(DIMV(rhoh))

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(rhoh)
      hi(1) = ARG_H1(rhoh)
      lo(2) = ARG_L2(rhoh)
      hi(2) = ARG_H2(rhoh)
      lo(3) = ARG_L3(rhoh)
      hi(3) = ARG_H3(rhoh)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      klo = max(lo(3),domlo(3))
      khi = min(hi(3),domhi(3))
      
      call filcc (rhoh,DIMS(rhoh),domlo,domhi,delta,xlo,bc)
      
c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(XHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(YLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(YHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(ZLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(ZHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "INFL_FORCE_F.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer DIMDEC(uflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T  uflct(:,:,:), vflct(:,:,:),  wflct(:,:,:)
      allocatable uflct, vflct, wflct

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(vel)
      hi(1) = ARG_H1(vel)
      lo(2) = ARG_L2(vel)
      hi(2) = ARG_H2(vel)
      lo(3) = ARG_L3(vel)
      hi(3) = ARG_H3(vel)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      klo = max(lo(3),domlo(3))
      khi = min(hi(3),domhi(3)) 

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

      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),1),
     $            DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,1))
      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),2),
     $            DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,2))
      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),3),
     $            DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,3))
      
c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.

      if (lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (bc(1,1,1).eq.EXT_DIR) then
                     if (forceLo .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,1) = u + uflct(1,j,k)
                     else
                        vel(i,j,k,1) = u
                     endif
                  endif
                  if (bc(1,1,2).eq.EXT_DIR) then
                     if (forceLo .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,2) = v + vflct(1,j,k)
                     else
                        vel(i,j,k,2) = v
                     endif
                  endif
                  if (bc(1,1,3).eq.EXT_DIR) then
                     if (forceLo .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,3) = w + wflct(1,j,k)
                     else
                        vel(i,j,k,3) = w
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (bc(1,2,1).eq.EXT_DIR) then
                     if (forceHi .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,1) = u + uflct(1,j,k)
                     else
                        vel(i,j,k,1) = u
                     endif
                  endif
                  if (bc(1,2,2).eq.EXT_DIR) then
                     if (forceHi .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,2) = v + vflct(1,j,k)
                     else
                        vel(i,j,k,2) = v
                     endif
                  endif
                  if (bc(1,2,3).eq.EXT_DIR) then
                     if (forceHi .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,3) = w + wflct(1,j,k)
                     else
                        vel(i,j,k,3) = w
                     endif
                  endif
               enddo
            enddo
         enddo
      endif    

      if (lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (bc(2,1,1).eq.EXT_DIR) then
                     if (forceLo .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,1) = u + uflct(i,1,k)
                     else
                        vel(i,j,k,1) = u
                     endif
                  endif
                  if (bc(2,1,2).eq.EXT_DIR) then
                     if (forceLo .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,2) = v + vflct(i,1,k)
                     else
                        vel(i,j,k,2) = v
                     endif
                  endif
                  if (bc(2,1,3).eq.EXT_DIR) then
                     if (forceLo .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,3) = w + wflct(i,1,k)
                     else
                        vel(i,j,k,3) = w
                     endif
                  endif
               enddo
            enddo
         enddo
      endif    
      
      if (hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(YHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (bc(2,2,1).eq.EXT_DIR) then
                     if (forceHi .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,1) = u + uflct(i,1,k)
                     else
                        vel(i,j,k,1) = u
                     endif
                  endif
                  if (bc(2,2,2).eq.EXT_DIR) then
                     if (forceHi .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,2) = v + vflct(i,1,k)
                     else
                        vel(i,j,k,2) = v
                     endif
                  endif
                  if (bc(2,2,3).eq.EXT_DIR) then
                     if (forceHi .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,3) = w + wflct(i,1,k)
                     else
                        vel(i,j,k,3) = w
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (bc(3,1,1).eq.EXT_DIR) then
                     if (forceLo .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,1) = u + uflct(i,j,1)
                     else
                        vel(i,j,k,1) = u
                     endif
                  endif
                  if (bc(3,1,2).eq.EXT_DIR) then
                     if (forceLo .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,2) = v + vflct(i,j,1)
                     else
                        vel(i,j,k,2) = v
                     endif
                  endif
                  if (bc(3,1,3).eq.EXT_DIR) then
                     if (forceLo .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,3) = w + wflct(i,j,1)
                     else
                        vel(i,j,k,3) = w
                     endif
                  endif
               enddo
            enddo
         enddo
      endif    
      
      if (hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(ZHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (bc(3,2,1).eq.EXT_DIR) then
                     if (forceHi .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,1) = u + uflct(i,j,1)
                     else
                        vel(i,j,k,1) = u
                     endif
                  endif
                  if (bc(3,2,2).eq.EXT_DIR) then
                     if (forceHi .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,2) = v + vflct(i,j,1)
                     else
                        vel(i,j,k,2) = v                        
                     endif
                  endif
                  if (bc(3,2,3).eq.EXT_DIR) then
                     if (forceHi .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,3) = w + wflct(i,j,1)
                     else
                        vel(i,j,k,3) = w
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (allocated(uflct)) then
         deallocate(uflct)
         deallocate(vflct)
         deallocate(wflct)
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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "INFL_FORCE_F.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer DIMDEC(uflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T uflct(:,:,:)
      allocatable uflct

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(xvel)
      hi(1) = ARG_H1(xvel)
      lo(2) = ARG_L2(xvel)
      hi(2) = ARG_H2(xvel)
      lo(3) = ARG_L3(xvel)
      hi(3) = ARG_H3(xvel)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      klo = max(lo(3),domlo(3))
      khi = min(hi(3),domhi(3)) 

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
      endif

      call filcc (xvel,DIMS(xvel),domlo,domhi,delta,xlo,bc)
      
c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceLo .and. strmwse_dir .eq. 1) then
                     xvel(i,j,k) = u + uflct(1,j,k)
                  else
                     xvel(i,j,k) = u
                  endif
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
                  call bcfunction(XHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceHi .and. strmwse_dir .eq. 1) then
                     xvel(i,j,k) = u + uflct(1,j,k)
                  else
                     xvel(i,j,k) = u
                  endif
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
                  call bcfunction(YLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceLo .and. strmwse_dir .eq. 2) then
                     xvel(i,j,k) = u + uflct(i,1,k)
                  else
                     xvel(i,j,k) = u
                  endif
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
                  call bcfunction(YHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceHi .and. strmwse_dir .eq. 2) then
                     xvel(i,j,k) = u + uflct(i,1,k)
                  else
                     xvel(i,j,k) = u
                  endif
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
                  call bcfunction(ZLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceLo .and. strmwse_dir .eq. 3) then
                     xvel(i,j,k) = u + uflct(i,j,1)
                  else
                     xvel(i,j,k) = u
                  endif
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
                  call bcfunction(ZHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceHi .and. strmwse_dir .eq. 3) then
                     xvel(i,j,k) = u + uflct(i,j,1)
                  else
                     xvel(i,j,k) = u
                  endif
               enddo
            enddo
         enddo
      endif

      if (allocated(uflct)) deallocate(uflct)

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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "INFL_FORCE_F.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
      integer DIMDEC(vflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T vflct(:,:,:)
      allocatable vflct

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(yvel)
      hi(1) = ARG_H1(yvel)
      lo(2) = ARG_L2(yvel)
      hi(2) = ARG_H2(yvel)
      lo(3) = ARG_L3(yvel)
      hi(3) = ARG_H3(yvel)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      klo = max(lo(3),domlo(3))
      khi = min(hi(3),domhi(3))
      
      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         enddo
         loFlctArray(strmwse_dir) = 1
         hiFlctArray(strmwse_dir) = 1
         call SET_ARGS(DIMS(vflct), loFlctArray, hiFlctArray)
         allocate(vflct(DIMV(vflct)))
         call INFL_FILL(FLCT_YVEL, DIMS(vflct), vflct, xlo, delta, time,
     $                  bc, domnlo, domnhi)
      endif

      call filcc (yvel,DIMS(yvel),domlo,domhi,delta,xlo,bc)
c      
c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.
c
      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceLo .and. strmwse_dir .eq. 1) then
                     yvel(i,j,k) = v + vflct(1,j,k)
                  else
                     yvel(i,j,k) = v
                  endif
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
                  call bcfunction(XHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceHi .and. strmwse_dir .eq. 1) then
                     yvel(i,j,k) = v + vflct(1,j,k)
                  else
                     yvel(i,j,k) = v
                  endif
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
                  call bcfunction(YLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceLo .and. strmwse_dir .eq. 2) then
                     yvel(i,j,k) = v + vflct(i,1,k)
                  else
                     yvel(i,j,k) = v
                  endif
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
                  call bcfunction(YHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceHi .and. strmwse_dir .eq. 2) then
                     yvel(i,j,k) = v + vflct(i,1,k)
                  else
                     yvel(i,j,k) = v
                  endif
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
                  call bcfunction(ZLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceLo .and. strmwse_dir .eq. 3) then
                     yvel(i,j,k) = v + vflct(i,j,1)
                  else
                     yvel(i,j,k) = v
                  endif
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
                  call bcfunction(ZHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceHi .and. strmwse_dir .eq. 3) then
                     yvel(i,j,k) = v + vflct(i,j,1)
                  else
                     yvel(i,j,k) = v
                  endif
               enddo
            enddo
         enddo
      endif

      if (allocated(vflct)) deallocate(vflct)

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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "INFL_FORCE_F.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
      integer DIMDEC(wflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T wflct(:,:,:)
      allocatable wflct

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(zvel)
      hi(1) = ARG_H1(zvel)
      lo(2) = ARG_L2(zvel)
      hi(2) = ARG_H2(zvel)
      lo(3) = ARG_L3(zvel)
      hi(3) = ARG_H3(zvel)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      klo = max(lo(3),domlo(3))
      khi = min(hi(3),domhi(3))


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
      
      call filcc (zvel,DIMS(zvel),domlo,domhi,delta,xlo,bc)
      
c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceLo .and. strmwse_dir .eq. 1) then
                     zvel(i,j,k) = w + wflct(1,j,k)
                  else
                     zvel(i,j,k) = w
                  endif
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
                  call bcfunction(XHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceHi .and. strmwse_dir .eq. 1) then
                     zvel(i,j,k) = w + wflct(1,j,k)
                  else
                     zvel(i,j,k) = w
                  endif
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
                  call bcfunction(YLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceLo .and. strmwse_dir .eq. 2) then
                     zvel(i,j,k) = w + wflct(i,1,k)
                  else
                     zvel(i,j,k) = w
                  endif
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
                  call bcfunction(YHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceHi .and. strmwse_dir .eq. 2) then
                     zvel(i,j,k) = w + wflct(i,1,k)
                  else
                     zvel(i,j,k) = w
                  endif
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
                  call bcfunction(ZLO, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceLo .and. strmwse_dir .eq. 3) then
                     zvel(i,j,k) = w + wflct(i,j,1)
                  else
                     zvel(i,j,k) = w
                  endif
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
                  call bcfunction(ZHI, x, y, z, time, u, v, w, rho, Yl, T, h)
                  if (forceHi .and. strmwse_dir .eq. 3) then
                     zvel(i,j,k) = w + wflct(i,j,1)
                  else
                     zvel(i,j,k) = w
                  endif
               enddo
            enddo
         enddo
      endif

      if (allocated(wflct)) deallocate(wflct)

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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(rhoY)
      hi(1) = ARG_H1(rhoY)
      lo(2) = ARG_L2(rhoY)
      hi(2) = ARG_H2(rhoY)
      lo(3) = ARG_L3(rhoY)
      hi(3) = ARG_H3(rhoY)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      klo = max(lo(3),domlo(3))
      khi = min(hi(3),domhi(3))
      
      call filcc (rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,bc)
      
c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(XLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(XHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(YLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(YHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(ZLO, x, y, z, time, u, v, w, rho, Yl, T, h)
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
                  call bcfunction(ZHI, x, y, z, time, u, v, w, rho, Yl, T, h)
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
      ihi = min(ARG_H1(p),domhi(1))
      jlo = max(ARG_L2(p),domlo(2))
      jhi = min(ARG_H2(p),domhi(2))
      klo = max(ARG_L3(p),domlo(3))
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

      end

      subroutine FORT_CALCDIVU(divu,DIMS(divu),plo,phi,xlo,dx,time)
c
c     This routine is for user defined divu, for debugging purposes only.
c
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
#include "cdwrk.H"
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
      
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               rad(i,j,k) = zero
            end do
         end do
      end do
      end
      
