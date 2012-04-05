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
#include "HEATTRANSFER_F.H"
#if defined(BL_DO_FLCT)
#include "infl_frc.H"
#endif

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

#define SDIM 3
#define USE_PMF_AVG

#if defined(BL_T3E) || defined(BL_AIX) || defined(BL_IRIX64) || defined(BL_CRAYX1)
#define OPEN_APPEND position='append'
#else
#define OPEN_APPEND access='append'
#endif

      subroutine FORT_ACTIVECONTROL(coft,time,dt,myproc,step,restart)

      implicit none

#include "probdata.H"
#include "cdwrk.H"
#include "bc.H"

      REAL_T coft,time,dt,vslope,slocal,V_new,dVmax
      integer myproc,step,ierr,restart
      REAL_T r1,r2,r3,r4,r5,r6,r7
      integer i1
      logical found_it

      if (restart.ne.0) then
         open(13,file=ac_hist_file,form='formatted',
     &        status='old',iostat=ierr)
         found_it = .false.
         if (ierr .eq. 0) then
            if (myproc.eq.0) then
               print*, 'Setting active control from history file ...'
            endif
            rewind(13)
            do
c
c                 This read(13) must correspond to the below write(13)
c
               read(13,1000,iostat=ierr) i1,r1,r2,r3,r4,r5,r6,r7
               if (ierr.ne.0) goto 100
               if (i1.eq.step) then

                  found_it = .true.
                  Vin = r2
                  tbase_control = r3
                  zbase_control = r4
                  dV_control = r5
                  sest = r6
                  coft_old = r7

               endif
            enddo

         else

            if (myproc.eq.0) then
               open(13,file=ac_hist_file,form='formatted', status='new')
            endif

         endif

 100     if (found_it .eqv. .false.) then

            if (myproc.eq.0) then
               print*, 'Setting active control to defaults ...'
            endif

         end if
         close(13)
         return
      end if

      if (coft_old.lt.zero) coft_old = coft

      if (dt.gt.zero) then
         zbase_control = zbase_control + Vin*dt + dV_control*dt**2
         V_in_old = Vin
         Vin = Vin + dt*dV_control

         slocal = half*(V_in_old + Vin) - (coft - coft_old)/(dt*scale_control)
         sest = (one - corr)*sest + corr*slocal
         vslope = two*((cfix-coft)/(scale_control*tau_control) + sest - Vin)/tau_control
         
         V_new = Vin + dt*vslope
         dVmax = changeMax_control * one
         V_new = MIN(MAX(V_new,Vin-dVmax),Vin+dVmax)
         V_new = MAX(zero,V_new)
         
         tbase_control = time
         dV_control = (V_new - Vin)/dt
         
         if (myproc.eq.0) then
            print *
            print *,'****************** control:'
            print *,'time,dt,coft,V_new:',time,dt,coft,V_new
            print *,'changeMax_control: ', changeMax_control
            print *,'vslope,sest,cfix,slocal:',vslope,sest,cfix,slocal
            print *,'coft_old,V_in_old:',coft_old,V_in_old
         endif
         
         coft_old = coft
         
         if (myproc.eq.0) then
            open(13,file=ac_hist_file,form='formatted',OPEN_APPEND)
            write(13,1000) step,time,Vin,tbase_control,zbase_control,
     &           dV_control,sest,coft_old
            close(13)
         endif
      endif
 1000 format(i7,7g26.18)
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
      call bl_abort('Should not be in FORT_HACK for this problem...')
      end

      subroutine FORT_SET_PROB_SPEC(fuel, oxid, prod, numspec)
      implicit none
#include "probdata.H"
      integer fuel, oxid, prod, numspec,isioproc
      call bl_pd_is_ioproc(isioproc)
      fuelID = fuel + 1
      oxidID = oxid + 1
      prodID = prod + 1
      if (isioproc.eq.1) then
         print*, 'fuelID: ', fuelID, ', oxidID: ', oxidID, 'prodID: ', prodID
      end if
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
#include "cdwrk.H"
#include "htdata.H"
#include "bc.H"
#if defined(BL_DO_FLCT)
#include "INFL_FORCE_F.H"
#endif
#include "visc.H"
#include "conp.H"

      integer lo(SDIM), hi(SDIM), i, j, zone, iCH4, getZone
      character*(maxspnml) spName
      INTEGER dimFile(3)
      REAL_T dxFile(3)
      REAL_T rand

      namelist /fortin/ vorterr, temperr, adverr, tempgrad, 
     &                  flametracval, probtype,
     &		        max_temp_lev, max_vort_lev, max_trac_lev,
     &                  traceSpecVal, nozzle_alt_lev,
     &                  max_nozzle_lev, refine_nozzle, refine_nozzle_x, 
     &                  refine_nozzle_y, refine_nozzle_z, refine_nozzle_z_alt, nozzle_width,
     &                  max_stick_lev, refine_stick, refine_stick_x, 
     &                  refine_stick_z, direction, min_nogrid_lev, nogrid_radius,
#if defined(BL_DO_FLCT)
     &                  forceInflow, numInflPlanesStore, strmwse_dir,
     $                  forceLo, forceHi, nCompInflow, flct_file, tstart_turb, 
#endif
     $                  turb_scale, Vin, Tin, phiin, rhot, zhot,
     $                  zstandoff, Thot, swK, swW, stTh, Ro, Rf, wallTh,
     &                  zBL, dBL, stBL, hotspot_r, hotspot_z, hotspot_dz,
     &                  hotspot_T, hotspot_htc, hotspot_w, vheight,
     &                  Vco_l,Vco_r,tVco_l,tVco_r,anisotsc,
     &                  InletFact_l,InletFact_r,tInletFact_l,tInletFact_r,
     &                  InletFluctFact_l,InletFluctFact_r,tInletFluctFact_l,tInletFluctFact_r,
     &                  tau_control, sest, cfix, changeMax_control, h_control,
     &			zbase_control, tbase_control,
     &                  thickFacTR, thickFacCH,
     &                  wire_x, wire_y, wire_z, wire_w, wire_h, wire_T, wire_r,
     &                  Tpipe, Tstick, stL, stLw, zInletEval, pertmag,
     &                  Rtran, Rftran, zfuelinit, Rfu, nchemdiag
      namelist /heattransin/ pamb, dpdt_factor
c
c      Build `probin' filename -- the name of file containing fortin namelist.
c
      integer maxlen, nCompFile, proc, isioproc
      parameter (maxlen=256)

      character probin*(maxlen)
      REAL_T y, yEval, pmf_vals(maxspec+3)
      integer n

      call bl_pd_myproc(proc)

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
      
      pertmag = zero
      zInletEval = zero
      anisotsc = 1.3d0
      vorterr = 1.e20
      temperr = zero
      adverr = 1.e20
      tempgrad  = 50.0d0
      flametracval = 0.0001d0
      probtype = 1
      max_temp_lev = 0
      max_vort_lev = 0
      max_trac_lev = 100
      traceSpecVal = 1.d-10
      refine_nozzle = 0
      refine_nozzle_z_alt = zero
      refine_nozzle_z = zero
      refine_nozzle_y = zero
      refine_nozzle_x = zero
      nozzle_width = -one
      max_nozzle_lev = 100
      nozzle_alt_lev = 100
      refine_stick = 0
      refine_stick_x = zero
      refine_stick_z = zero
      max_stick_lev = 100
      pamb = 101325.d0
      dpdt_factor = 0.3d0
      direction = "z"
      vheight = zero
      min_nogrid_lev = 100
      nogrid_radius = 1.e10

      add_turb = .FALSE.

#if defined(BL_DO_FLCT)
      forceInflow = .FALSE.
      numInflPlanesStore = -1
      numInflPlanesStore = 30
c
c     Don't need to default 'nCompInflow' as it is block data'd to /3/
c
      forceLo = .TRUE.
      forceHi = .FALSE.
      strmwse_dir = FLCT_ZVEL
      flct_file = ''
#endif
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
      Tin = -1.d0
      phiin = -1.d0

      Vco_l  = 0.1d0
      Vco_r  = 0.1d0
      tVco_l = 0.0d0
      tVco_r = 0.0d0

      InletFact_l = -1
      InletFact_r = -1
      tInletFact_l = -1
      tInletFact_r = -1

      InletFluctFact_l = -1
      InletFluctFact_r = -1
      tInletFluctFact_l = -1
      tInletFluctFact_r = -1

      Thot = 2.d3
      rhot = Rf
      zhot = Rf
      zstandoff = zero
      swK = zero
      swW = 0.6d0

c     Initialize control variables
      tau_control = one
      sest = zero
      corr = one
      changeMax_control = .05
      coft_old = -one
      cfix = zero
      ac_hist_file = 'AC_History'
      dV_control = zero
      tbase_control = zero
      h_control = -1

      wire_x = -1.0
      wire_y = -1.0
      wire_z = -1.0
      wire_w = -1.0
      wire_h = zero
      wire_T = -1.0
      wire_r = -1.0
      zbase_control = zero

c     If the following two Temperatures do not appear in probin, they default to the 
c     compiled-in values in the data statements in setupbc
      Tpipe = -1
      Tstick = -1

c     Length of stick at Tstick, then ramps to Tpipe over layer stLw
      stL = 20*Ro
      stLw = .1*Ro

c     New swirl parameters
      Rtran = .1d0
      Rftran = .1d0
      zfuelinit = .1d0
      Rfu = .1d0
      
      nchemdiag = 1

      read(untin,fortin)

c     Initialize control variables that depend on fortin variables
      V_in_old = Vin

      if (max_vort_lev.lt.0) max_vort_lev=max_temp_lev
      
      read(untin,heattransin)

   10 close(unit=untin)

#if defined(BL_DO_FLCT)
      if (forceInflow .eqv. .FALSE.) then
         forceLo = .FALSE.
         forceHi = .FALSE.
      else
         if (flct_file.ne.'') then
#define FF_UNIT 20
            if (proc.eq.0) print*, '...initializing turbulence'
            open(FF_UNIT, file=trim(flct_file)//'/HDR',form='formatted',status='old')
            call RD_SCL_FLCTHD(FF_UNIT,nCompFile,dimFile,probSizeFile,dxFile)
            close(FF_UNIT)
         endif
      endif

c     We build our our "time" value in VELFILL for stepping through the
c     turbulence file.  The routine relies on time not being further scaled
c     by convVel within the INFL stuff....set convVel = 1.0     
      convVel = one
#endif

c     Initialize the control parameters...even if we are not controlling, as they are
c     used in the calculation of the time offset into the turbulence file
      dV_control = zero

c     Load domain dimensions into common, and set up boundary functions
      domnlo(1) = problo(1)
      domnlo(2) = problo(2)
      domnlo(3) = problo(3)
      domnhi(1) = probhi(1)
      domnhi(2) = probhi(2)
      domnhi(3) = probhi(3)
      
      call setupbc()
      bcinit = .true.
      
      iCH4 = -1
      do n = 1,Nspec
         call get_spec_name(spName,n)
         if (spName .eq. 'CH4') iCH4 = n-1
      end do
      zone = getZone(half*(domnlo(1)+domnhi(1)),half*(domnlo(2)+domnhi(2)),domnlo(3))
      scale_control = Y_bc(iCH4,zone)*rho_bc(zone)*(domnhi(1)-domnlo(1))*(domnhi(2)-domnlo(2))
      if (h_control .gt. zero) then
         cfix = scale_control * h_control
      endif

      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
      end if

      end

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
#define ZLO 4
#define ZHI 5

      subroutine setupbc()
      implicit none
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "htdata.H"
      
      REAL_T Patm, x, pmf_vals(maxspec+3), sum, Y_pmf(maxspec), a
      integer i, zone, n, lo(SDIM), hi(SDIM), iN2, iO2, getZone
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

      integer isioproc


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
c     phi = 0.8
c     data X_CH4 /  0.0774907749077d0,   zero, 0.0774907749077d0,   zero,    zero,   zero ,   zero /
c     data X_O2  /   0.193726937269d0,  0.21d0, 0.193726937269d0,  0.21d0,  0.21d0, 0.21d0,  0.21d0 /
      data T     /        3.d2,          3.d2,         3.d2,        3.d2,    3.d2,   3.d2,   3.d2  /
#endif
      Patm = pamb / 101325.0d0
      iN2 = -1
      do n = 1,Nspec
         call get_spec_name(name, n)
         if (name .eq. 'N2' ) iN2 = n
      end do

      if (Tpipe.gt.0) T(BL_PIPEEND) = Tpipe
      if (Tstick.gt.0) T(BL_STICK) = Tstick

      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         print *, '*********** zone boundary values at startup'
      endif
      do zone = 1,Nzones
c         if(isioproc.eq.1)print *,"in zone loop",zone, BL_FUELPIPE
         if (zone.eq.BL_FUELPIPE .and. ((probtype.eq.4).or.(probtype.eq.5)) ) then
c         if(isioproc.eq.1)print *,"inside a",zone, BL_FUELPIPE
            if (phiin .ge. 0.d0) then
               iN2  = -1
               iO2  = -1
               do n = 1,Nspec
                  call get_spec_name(name,n)
                  if (name .eq. 'N2' ) iN2  = n
                  if (name .eq. 'O2' ) iO2  = n
               enddo
 
               do n=1,Nspec
                  Xt(n) = zero
               enddo

               call get_spec_name(name,fuelID)

               if (name .eq. 'CH4') then
c
c                 hc + a.O2 -> b.CO2 + c.H2O
c                 for hc = CH4, a=2, b=1, c=2
c
                  a = 2.0d0

                  Xt(iO2)    = 1.d0/(1.d0 + phiin/a  + 0.79d0/0.21d0)
                  Xt(fuelID) = phiin * Xt(iO2) / a
                  Xt(iN2)    = 1.d0 - Xt(fuelID) - Xt(iO2)
               else if (name .eq. 'H2') then
c
c                 hc + a.O2 -> b.CO2 + c.H2O
c                 for hc = H2, a=.5, b=0, c=1
c
                  a = .5d0

                  Xt(iO2)    = 1.d0/(1.d0 + phiin/a  + 0.79d0/0.21d0)
                  Xt(fuelID) = phiin * Xt(iO2) / a
                  Xt(iN2)    = 1.d0 - Xt(fuelID) - Xt(iO2)
               else if (name .eq. 'C3H8') then
c
c                 hc + a.O2 -> b.CO2 + c.H2O
c                 for hc = C3H8, a=5, b=3, c=4
c
                  a = 5.0d0

                  Xt(iO2)    = 1.d0/(1.d0 + phiin/a  + 0.79d0/0.21d0)
                  Xt(fuelID) = phiin * Xt(iO2) / a
                  Xt(iN2)    = 1.d0 - Xt(fuelID) - Xt(iO2)
               else
                  if (isioproc.eq.1) then
                     print*, 'setupbc(): unknown fuel type'
                  end if
                  call bl_abort('')
               end if

               T_bc(zone) = Tin
            else
c            if(isioproc.eq.1)print *,"inside b",zone, BL_FUELPIPE
               call pmf(zInletEval,zInletEval,pmf_vals,n)
               if (n.ne.Nspec+3) then
                  call bl_abort('setupbc: n .ne. Nspec+3')
               endif
               do n =1,Nspec
                  Xt(n) = pmf_vals(3+n)
               enddo
               T_bc(zone) = pmf_vals(1)
            endif
         else
c         if(isioproc.eq.1)print *,"inside c",zone, BL_FUELPIPE
            sum = zero
            do n = 1,Nspec
               Xt(n) = zero
               call get_spec_name(name, n)
               if (name .EQ. 'CH4') Xt(n) = X_CH4(zone)
               if (name .EQ. 'O2')  Xt(n) = X_O2(zone)
               if (name .NE. 'N2')  sum = sum + Xt(n)
            end do
            Xt(iN2) = one - sum
            T_bc(zone) = T(zone)
         endif
         CALL CKXTY (Xt, IWRK, RWRK, Y_bc(0,zone))

         call FORT_RHOfromPTY(lo,hi,
     &        rho_bc(zone),ARLIM(lo),ARLIM(hi),
     &        T_bc(zone),ARLIM(lo),ARLIM(hi),
     &        Y_bc(0,zone),ARLIM(lo),ARLIM(hi), Patm)
         call FORT_HMIXfromTY(lo,hi,
     &        h_bc(zone),ARLIM(lo),ARLIM(hi),
     &        T_bc(zone),ARLIM(lo),ARLIM(hi),
     &        Y_bc(0,zone),ARLIM(lo),ARLIM(hi))

         if (isioproc.eq.1) then
            print *, 'zone',zone
            print *, '  T,rho,h: ',T_bc(zone),rho_bc(zone),h_bc(zone)
            print *, '  Y:'
            print *, (Y_bc(i,zone),i=0,Nspec-1)
         endif

      end do

      end

      integer function getZone(x, y, z)
      implicit none
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      REAL_T x, y, z, r, x1
c
c     ZONE DEF:
c     BL_FUELPIPE, BL_COFLOW, BL_STICK, BL_WALL, BL_AMBIENT, BL_VOLUME
c
      getZone = BL_VOLUME
      if ( ((direction.eq."z") .and. (z.le.domnlo(3))).or.
     &     ((direction.eq."y") .and. (y.le.domnlo(2))).or.
     &     ((direction.eq."x") .and. (x.le.domnlo(1)))) then

         if ( (probtype.ge.4) .or. (probtype.eq.5) ) then

            getZone = BL_FUELPIPE

         else

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
     &     Vf,Vco_l,Vco_r,tVco_l,tVco_r,time,swK,swW,probtype,coCUTOFFp)
      implicit none
      REAL_T V(SDIM)
      REAL_T x,y,dx(SDIM),Rf,wallTh,Ro,stTh,dBL,stBL
      REAL_T Vf,Vco_l,Vco_r,tVco_l,tVco_r,time,swK,swW,coCUTOFFp,coCUTOFF
      integer probtype

      REAL_T incx,incy,x1,x2,eta
      REAL_T Vincx,Vincy,Vincz,Vswirl,r
      integer i,j,M
      parameter (M=10)
c
c     size of subvolumes, there are M*M of them
c
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
               Vincz = Vf*TANH(two*(one - r/Ro)/dBL)
#endif
            else if (r.ge.Ro+wallTh) then
               if (time .le. tVco_l) then
                  Vincz = Vco_l
               else if (time .ge. tVco_r) then
                  Vincz = Vco_r
               else
                  Vincz = Vco_l+(time-tVco_l)*(Vco_r-Vco_l)/(tVco_r-tVco_l)
               endif

c     For stagnation case, ramp co-flow to zero at coCUTOFF
c              if (probtype.eq.3) then
c                 coCUTOFF = coCUTOFFp - Ro - wallTh
c                 Vincz = Vincz*MAX(zero,(coCUTOFF-r+Ro+wallTh)/coCUTOFF)
c              endif

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
c              Vswirl = Vf*(-0.067641d0 + 2.0113d0*eta - 19.305d0*eta**2 + 72.404d0*eta**3
c     &                     - 126.14d0*eta**4 + 111.66d0*eta**5 - 31.098d0*eta**6)
c              Vradial =Vf*(-0.0022569d0 + 0.3562d0*eta - 0.84751d0*eta**2 + 1.8533d0*eta**3
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
      REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h, dx(SDIM), bin, bout, Tout
      REAL_T rr
      logical getuvw

      integer isioproc

#include "cdwrk.H"
#include "htdata.H"
#include "bc.H"
#include "probdata.H"

      REAL_T Vbndry(3), eta, Patm, r, rn, u_inf(3), Vco
      integer n, getZone, zone, lo(SDIM), hi(SDIM)
      data lo / 1, 1, 1 /
      data hi / 1, 1, 1 /

      double precision factor

      call bl_pd_is_ioproc(isioproc)

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if
      
      if (orient .le. ZHI) then
         
         if (probtype.eq.5) then
            
            r = SQRT( x*x + y*y )
            
c     Rftran is the transition width between fuel and air
            eta =  0.5d0*(1.d0 - TANH(2*(r-Rfu)/Rftran))
     &           * 0.5d0*(1.d0 - TANH(2*(z-zfuelinit)/Rftran))
            
            do n = 0, Nspec-1
               Yl(n) = (1.d0-eta)*Y_bc(n,BL_COFLOW) + eta*Y_bc(n,BL_FUELPIPE)
            end do
            T = (1.d0-eta)*T_bc(BL_COFLOW) + eta*T_bc(BL_FUELPIPE)

            rr = SQRT( x*x + y*y + (z-zhot)**2 )
            eta = 0.5d0*(1.d0 - TANH(2*(rr-rhot)/Rftran))
            do n = 0, Nspec-1
               Yl(n) = (1.d0-eta)*Yl(n) + eta*Y_bc(n,BL_COFLOW)
            end do
            T = (1.d0-eta)*T + eta*Thot

            Patm = pamb / 101325.0d0
            
            call FORT_RHOfromPTY(lo,hi,
     &           rho,ARLIM(lo),ARLIM(hi),
     &           T,ARLIM(lo),ARLIM(hi),
     &           Yl,ARLIM(lo),ARLIM(hi), Patm)
            call FORT_HMIXfromTY(lo,hi,
     &           h,ARLIM(lo),ARLIM(hi),
     &           T,ARLIM(lo),ARLIM(hi),
     &           Yl,ARLIM(lo),ARLIM(hi))

            if (getuvw) then

c     Compute coflow velocity components (allows ramping up/down)
               u_inf(1) = 0.d0
               u_inf(2) = 0.d0
               if ( (time .gt. zero).and.(tVco_l.ge.zero) ) then
                  if (time .le. tVco_l) then
                     u_inf(3) = Vco_l
                  else if (time .ge. tVco_r) then
                     u_inf(3) = Vco_r
                  else
                     u_inf(3) = Vco_l+(time-tVco_l)*(Vco_r-Vco_l)/(tVco_r-tVco_l)
                  endif
               endif

c     Load experimental database values
#ifdef SWIRL_DATA
               call vswirlXYZ(x,y,u,v,w)
#else
               u = 0
               v = 0
               w = 0
#endif
               
c     Compute a time-varying factor to apply to the velocity fields above the fuel pipe
               factor = 1.d0
               if ( (time .gt. zero).and.(tInletFact_l.ge.zero) ) then
                  if (time .le. tInletFact_l) then
                     factor = InletFact_l
                  else if (time .ge. tInletFact_r) then
                     factor = InletFact_r
                  else
                     factor = InletFact_l
     &                    +(time-tInletFact_l)*(InletFact_r-InletFact_l)
     &                    / (tInletFact_r-tInletFact_l)
                  endif
               endif

               u = u*factor
               v = v*factor
               w = w*factor

c     Continue expt data with TANH profiles to u_inf values
               if (r.gt.Rfu) then
                  eta = TANH(2*(r-Rfu)/Rtran)
                  u = eta*u_inf(1) + (1.d0-eta)*u
                  v = eta*u_inf(2) + (1.d0-eta)*v
                  w = eta*u_inf(3) + (1.d0-eta)*w
               endif
            endif

         else

            zone = getZone(x,y,z)
            rho = rho_bc(zone)
            do n = 0, Nspec-1
               Yl(n) = Y_bc(n,zone)
            end do
            T = T_bc(zone)
            h = h_bc(zone)
            Patm = pamb / 101325.0d0
            
            r = SQRT(x**2 + y**2)
            
            if (probtype.ne.4 .and.
     &           zone.eq.BL_FUELPIPE  .or.
     &           zone.eq.BL_PIPEEND   .or.
     &           zone.eq.BL_COFLOW    .or.
     &           zone.eq.BL_STICK          ) then

               T = T_bc(BL_FUELPIPE)

               if (stTh.gt.zero) then
                  
c     bin is the blend from the stick center T to T of the pipe end over stLw
                  bin = one - half*(one+TANH((ABS(y)-half*stL)/(fourth*stLw)))
                  T = T_bc(BL_STICK)*bin + T_bc(BL_PIPEEND)*(one-bin)
                  
c     now, bin is the blend from the stick to the fuel or pipe end over stBL
                  if (zone.eq.BL_FUELPIPE.or.zone.eq.BL_PIPEEND.or.zone.eq.BL_STICK) then
                     if (r.le.Ro) then
                        Tout = T_bc(BL_FUELPIPE)
                     else
                        Tout = T_bc(BL_PIPEEND)
                     endif
                     bin = one - half*(one+TANH((ABS(x)-half*stTh)/(fourth*stBL)))
                     T = T*bin + Tout*(one - bin)
                  endif
                  
               endif
               
c     now, bin is the blend from fuel stream to pipe end, bout is the blend from pipe end to coflow
               bin = one - half*(one+TANH((r-Ro)/(fourth*.3*wallTh)))
               bout = one - half*(one+TANH((r-(Ro+wallTh))/(fourth*.3*wallTh)))
               
               T = (T*bin + T_bc(BL_PIPEEND)*(one-bin))*bout + T_bc(BL_COFLOW)*(one-bout)
               
               if (zone.eq.BL_FUELPIPE) then
#ifdef SWIRL
                  rn = r/Ro
                  
c     c     Fit as of 4/29/03
c     eta = 9.920628558770E-001
c     &           -1.065813294748E-001 *rn**1
c     &           +1.800416951718E+000 *rn**2
c     &           -1.158582069687E+001 *rn**3
c     &           +3.588317172715E+001 *rn**4
c     &           -5.816747423739E+001 *rn**5
c     &           +4.632391073835E+001 *rn**6
c     &           -1.450956235808E+001 *rn**7
c     c     Ramp down the outter value to zero, 6/2/03
c     eta = eta*0.5*(1-TANH((rn-0.92)/.03))
                  
c     eta = 1.002d0 - 0.098868d0*r + 0.84254d0*r**2 - 2.6231d0*r**3 +
c     $           3.0126d0*r**4 - 1.9876d0*r**5 + 0.52246d0*r**6
c     
                  
c      Fit by Tecplot for 288x288 swirl inlet data, 1/30/04 (MSD)
                  
                  eta = 1.004677322175D+000
     &                - 9.017625583830D-002*rn
     &                + 3.611171220764D-001*rn**2
     &                - 5.644888675141D-001*rn**3
               
                  eta = min(max(eta,0.d0),1.d0)
#else
                  eta = one
#endif
                  do n = 0, Nspec-1
                     Yl(n) = (1.d0-eta)*Y_bc(n,BL_COFLOW) + eta*Y_bc(n,BL_FUELPIPE)
                  end do
               endif
               
               call FORT_RHOfromPTY(lo,hi,
     &                              rho,ARLIM(lo),ARLIM(hi),
     &                              T,ARLIM(lo),ARLIM(hi),
     &                              Yl,ARLIM(lo),ARLIM(hi), Patm)
               call FORT_HMIXfromTY(lo,hi,
     &                              h,ARLIM(lo),ARLIM(hi),
     &                              T,ARLIM(lo),ARLIM(hi),
     &                              Yl,ARLIM(lo),ARLIM(hi))
            endif

            if (getuvw) then
               if (probtype.eq.4) then
                  u = zero
                  v = zero
                  w =  Vin + (time-tbase_control)*dV_control
               else
                  call VAvg(Vbndry,x,y,dx,Rf,wallTh,Ro,stTh,dBL,stBL,
     &                      Vin,Vco_l,Vco_r,tVco_l,tVco_r,time,swK,swW,probtype,domnhi(1))
                  u = Vbndry(1)
                  v = Vbndry(2)
                  w = Vbndry(3)
               endif
            endif
            
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
#include "probdata.H"
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
#include "probdata.H"
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

#include "cdwrk.H"
#include "htdata.H"
#include "bc.H"
#include "probdata.H"
#if defined(BL_DO_FLCT)
#include "INFL_FORCE_F.H"
#endif

      integer i, j, k, n, iN2, iO2
      REAL_T x, y, z, zt, u, v, w, rho, Yl(maxspec), Xl(maxspec), T, h, Patm
      REAL_T turbSclX, turbSclY, turbSclZ, convVelSAVE
      REAL_T eta, eta0, pmf_vals(maxspec+3), bg_vals(maxspec+3), r, y1, y2, theta
      REAL_T zpert,Lx,Ly
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
c         
c     Make a platter of data in XY, will march over z
c
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
         convVelSAVE = convVel
         convVel = 15.d0
         if (.false. .and. flct_file.ne.'') then      
c         if (flct_file.ne.'') then      
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
         convVel = convVelSAVE
#endif
      enddo            

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
                        call bl_abort('INITDATA: n .ne. Nspec+3')
                     endif
                  endif

                  eta = one
                  if (rhot .gt. zero) eta = half*(one + TANH(two*(r-rhot)/zBL))

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
                  else if (Nspec.eq.3) then
                     do n = 1,Nspec+3
                        pmf_vals(n) = bg_vals(n)
                     enddo
                  else
                     call pmf(y1,y1,pmf_vals,n)
                     if (n.ne.Nspec+3) then
                        call bl_abort('INITDATA: n .ne. Nspec+3')
                     endif
                  endif

                  eta = zero
                  if (rhot .gt. zero) eta = half*(one + TANH(two*(r-rhot)/zBL))

                  scal(i,j,k,Temp) = (one - eta)*pmf_vals(1) + eta*bg_vals(1)
                  do n = 1,Nspec
                     Xl(n) = (one - eta)*pmf_vals(3+n) + eta*bg_vals(3+n)
                  end do 

                  eta = zero
                  if (zhot .gt. zero) eta = half*(one + TANH(two*(z-zhot)/zBL))

                  scal(i,j,k,Temp) = (one - eta)*scal(i,j,k,Temp) + eta*bg_vals(1)
                  do n = 1,Nspec
                     Xl(n) = (one - eta)*Xl(n) + eta*bg_vals(3+n)
                  end do 

                  CALL CKXTY (Xl, IWRK, RWRK, Yl)
                  do n = 1,Nspec
                     scal(i,j,k,FirstSpec+n-1) = Yl(n)
                  end do
                  scal(i,j,k,Trac) = zero

                  zt = domnlo(3) - half*delta(3)

c     Propagate mean inflow profiles through domain
                  call bcfunction(XLO,x,y,zt,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  vel(i,j,k,1) = u
                  vel(i,j,k,2) = v
                  vel(i,j,k,3) = w                  
c                  if (z-domnlo(3) .LT. vheight) then
c                     vel(i,j,k,1) = zero
c                     vel(i,j,k,2) = zero
c                     vel(i,j,k,3) = zero
c                  end if
               end do
            end do
         end do
      else if (probtype .eq. 4) then
         do k= lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)            
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)

                  zpert = 0.d0
                  if (pertmag .gt. 0.d0) then
                     Lx = domnhi(1) - domnlo(1)
                     Ly = domnhi(2) - domnlo(2)
                     zpert = pertmag*(1.000 * sin(2*Pi*4*x/Lx)             * sin(2*Pi*5*y/Ly)
     &                              + 1.023 * sin(2*Pi*2*(x-.004598)/Lx)   * sin(2*Pi*4*(y-.0053765)/Ly)
     &                              + 0.945 * sin(2*Pi*3*(x-.00712435)/Lx) * sin(2*Pi*3*(y-.02137)/Ly) 
     &                              + 1.017 * sin(2*Pi*5*(x-.0033)/Lx)     * sin(2*Pi*6*(y-.018)/Ly) 
     &                                          + .982 * sin(2*Pi*5*(x-.014234)/Lx) )
                  endif

                  y1 = (z-half*delta(3)-zstandoff+zpert)*100.d0
                  y2 = (z+half*delta(3)-zstandoff+zpert)*100.d0

                  call pmf(y1,y2,pmf_vals,n)

                  if (n.ne.Nspec+3) then
                     call bl_abort('INITDATA: n .ne. Nspec+3')
                  endif

                  scal(i,j,k,Temp) = pmf_vals(1)
                  do n = 1,Nspec
                     Xl(n) = pmf_vals(3+n)
                  end do 

                  CALL CKXTY (Xl, IWRK, RWRK, Yl)

                  do n = 1,Nspec
                     scal(i,j,k,FirstSpec+n-1) = Yl(n)
                  end do
c                  scal(i,j,k,Trac) = zero
                  scal(i,j,k,Trac) = -z

               end do
            end do
         end do
      else if (probtype .eq. 5) then
         do k= lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)            
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)

                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  scal(i,j,k,Temp) = T
                  do n = 1,Nspec
                     scal(i,j,k,FirstSpec+n-1) = Yl(n)
                  end do
                  scal(i,j,k,Trac) = zero

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
      !
      ! This is a version of INITDATA in which the velocities,
      ! the temperature & the mass fractions are passed in.
      !
      subroutine FORT_INITDATANEWMECH(level,time,lo,hi,nscal,
     &     vel,scal,DIMS(state),press,DIMS(press),
     &     delta,xlo,xhi)
      implicit none
      integer  level, nscal, nvel
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(state)
      integer  DIMDEC(press)
      REAL_T   xlo(SDIM), xhi(SDIM)
      REAL_T   time, delta(SDIM)
      REAL_T   vel(DIMV(state),SDIM)
      REAL_T   scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))

#include "cdwrk.H"
#include "htdata.H"
#include "bc.H"
#include "probdata.H"

      integer i, j, k, n
      REAL_T Patm

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               scal(i,j,k,Trac) = zero
            end do
         end do
      end do

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
#include "cdwrk.H"
#include "htdata.H"
      integer i, j, k, n, Tid, RHid, YSid, YEid, ys, ye
      integer getZone
      logical do_T, do_RH, do_Y
      REAL_T xl, xr, xh, yb, yt, yh, z

      if (probtype.lt.4) then
         Tid  = Temp      - id + SDIM
         RHid = RhoH      - id + SDIM
         YSid = FirstSpec - id + SDIM
         YEid = LastSpec  - id + SDIM
         
         do_T  = (Tid  .GE. 1) .AND. (Tid  .LE. ncomp)
         do_RH = (RHid .GE. 1) .AND. (RHid .LE. ncomp)
         ys = MAX(YSid,1)
         ye = MIN(YEid,ncomp)
         do_Y = (ye - ys + 1) .GE. 1
c     
c     Do species, Temp, rhoH
c     
         if ((idir.EQ.2) .AND. (lo(3) .LE. domlo(3))
     &           .AND. (do_T .OR. do_RH .OR. do_Y) ) then
               
            k = lo(3)
            z = float(k)*dx(3)+domnlo(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  
                  xl = float(i)*dx(1)+domnlo(1) 
                  xr = (float(i)+1.d0)*dx(1)+domnlo(1) 
                  xh = 0.5d0*(xl+xr)
                  yb = float(j)*dx(2)+domnlo(2) 
                  yt = (float(j)+1.d0)*dx(2)+domnlo(2) 
                  yh = 0.5d0*(yb+yt)
                  
                  if ( (getZone(xl,yb,z).eq.BL_STICK) .OR.
     &                 (getZone(xh,yb,z).eq.BL_STICK) .OR.
     &                 (getZone(xr,yb,z).eq.BL_STICK) .OR.
     &                 (getZone(xl,yh,z).eq.BL_STICK) .OR.
     &                 (getZone(xh,yh,z).eq.BL_STICK) .OR.
     &                 (getZone(xr,yh,z).eq.BL_STICK) .OR.
     &                 (getZone(xl,yt,z).eq.BL_STICK) .OR.
     &                 (getZone(xh,yt,z).eq.BL_STICK) .OR.
     &                 (getZone(xr,yt,z).eq.BL_STICK) .OR.
     &                 (getZone(xl,yb,z).eq.BL_PIPEEND) .OR.
     &                 (getZone(xh,yb,z).eq.BL_PIPEEND) .OR.
     &                 (getZone(xr,yb,z).eq.BL_PIPEEND) .OR.
     &                 (getZone(xl,yh,z).eq.BL_PIPEEND) .OR.
     &                 (getZone(xh,yh,z).eq.BL_PIPEEND) .OR.
     &                 (getZone(xr,yh,z).eq.BL_PIPEEND) .OR.
     &                 (getZone(xl,yt,z).eq.BL_PIPEEND) .OR.
     &                 (getZone(xh,yt,z).eq.BL_PIPEEND) .OR.
     &                 (getZone(xr,yt,z).eq.BL_PIPEEND) ) then
                     
c                    if (do_T)  diff(i,j,k,Tid ) = 0.d0
c                    if (do_RH) diff(i,j,k,RHid) = 0.d0
                     if (do_Y) then
                        do n=ys,ye
                           diff(i,j,k,n) = 0.d0
                        enddo
                     endif
                     
                  endif
               end do
            end do
         endif
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
      REAL_T    x, y, z, r
      logical   in_refine_zone

#include "probdata.H"

      if (level.lt.max_trac_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),
     &                 ftrac(i,j,k,1).gt.flametracval)
               enddo
            enddo
         enddo
      endif

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
                  tag(i,j,k) = merge(set,tag(i,j,k),in_refine_zone)
               end do
            end do
         end do
      end if
      if ((level.lt.nozzle_alt_lev).and.(refine_nozzle.eq.1)) then
         do k = lo(3), hi(3)
            z = (float(k)+.5)*dx(3)+problo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*dx(2)+problo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*dx(1)+problo(1)
                  r = SQRT(x*x+y*y)
                  in_refine_zone = (r .le. refine_nozzle_x)
                  in_refine_zone = in_refine_zone .and.
     &                  (z - problo(3)) .le. refine_nozzle_z_alt
                  tag(i,j,k) = merge(set,tag(i,j,k),in_refine_zone)
               end do
            end do
         end do
      end if
      if ((probtype.lt.4).and.(level.lt.max_stick_lev).and.(refine_stick.eq.1)) then
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
                  tag(i,j,k) = merge(set,tag(i,j,k),in_refine_zone)
               end do
            end do
         end do
      end if
      if (level.ge.min_nogrid_lev-1) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*dx(2)+problo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*dx(1)+problo(1)
                  r = SQRT(x*x+y*y)
                  in_refine_zone = (r.gt.nogrid_radius)
                  tag(i,j,k) = merge(clear,tag(i,j,k),in_refine_zone)
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
                  ax = abs(temperature(i+1,j,k,1) - temperature(i,j,k,1))
                  ay = abs(temperature(i,j+1,k,1) - temperature(i,j,k,1))
                  az = abs(temperature(i,j,k+1,1) - temperature(i,j,k,1))
                  ax = MAX(ax,abs(temperature(i,j,k,1) - temperature(i-1,j,k,1)))
                  ay = MAX(ay,abs(temperature(i,j,k,1) - temperature(i,j-1,k,1)))
                  az = MAX(az,abs(temperature(i,j,k,1) - temperature(i,j,k-1,1)))
                  aerr = max(ax,ay,az)
                  tag(i,j,k) = merge(set,tag(i,j,k),aerr.ge.tempgrad)
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
                  tag(i,j,k) = merge(set,tag(i,j,k),
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

      integer    i,j,k
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(adv)
      lo(2) = ARG_L2(adv)
      lo(3) = ARG_L3(adv)
      hi(1) = ARG_H1(adv)
      hi(2) = ARG_H2(adv)
      hi(3) = ARG_H3(adv)

      call filcc (adv,DIMS(adv),domlo,domhi,delta,xlo,bc)

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            do j = lo(2),hi(2)
               do i = lo(1), hi(1)
                  adv(i,j,k) = time
               enddo
            enddo
         enddo
      endif    
      
c      call fillWithZeros(adv,DIMS(adv),domlo,domhi,bc)

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

#include "cdwrk.H"
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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#if defined(BL_DO_FLCT)
#include "INFL_FORCE_F.H"
#endif
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x, eta
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
      REAL_T  cs,sg,r,scaler,scalet
      REAL_T  turbSclX,turbSclY,turbSclZ,turbSclR,turbSclT,factor
      REAL_T frcfct

      REAL_T, parameter :: T0 = .1436577d0

#if defined(BL_DO_FLCT)
      integer DIMDEC(uflct)
      REAL_T t_flct, dt_flct
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

      frcfct = 1.0d0

c      if (probtype.eq.5) then
c         frcfct = 1.0d0 + .2d0 * sin(2*Pi*347*(time-T0))
c      endif

c     Compute a time-varying factor to apply to the velocity field fluctuations
      factor = 1.d0
      if ( (time .gt. zero).and.(tInletFluctFact_l.ge.zero) ) then
         if (time .le. tInletFluctFact_l) then
            factor = InletFluctFact_l
         else if (time .ge. tInletFluctFact_r) then
            factor = InletFluctFact_r
         else
            factor = InletFluctFact_l
     &           +(time-tInletFluctFact_l)*(InletFluctFact_r-InletFluctFact_l)
     &           /(tInletFluctFact_r-tInletFluctFact_l)
         endif
      endif

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
c     Note that we are 'scaling time' here to step into the fluct file to the
c     correct depth.  This requires that time is not further scaled inside the
c     the INFL_FILL routine.  Just to be sure, we set convVel = 1 here again.
         convVel = one
         dt_flct = time - tbase_control
        t_flct = zbase_control + Vin*dt_flct + dV_control*dt_flct**2
c         t_flct = zbase_control + Vin*dt_flct*factor + dV_control*dt_flct**2

         call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo, delta, t_flct,
     $                  bc(1,1,1), domnlo, domnhi)
         call INFL_FILL(FLCT_YVEL, DIMS(uflct), vflct, xlo, delta, t_flct,
     $                  bc(1,1,2), domnlo, domnhi)
         call INFL_FILL(FLCT_ZVEL, DIMS(uflct), wflct, xlo, delta, t_flct,
     $                  bc(1,1,3), domnlo, domnhi)

      endif
#endif

      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),1),
     $            DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,1))
      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),2),
     $            DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,2))
      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),3),
     $            DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,3))

      if (lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5d0)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5d0)*delta(2)+domnlo(2)

                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(1,1,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,1) = u + uflct(1,j,k)*turbSclX(x,y,z)*turb_scale
                     endif
#endif
                  endif

                  if (bc(1,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,2) = v + vflct(1,j,k)*turbSclY(x,y,z)*turb_scale
                     endif
#endif
                  endif

                  if (bc(1,1,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,3) = w + wflct(1,j,k)*turbSclZ(x,y,z)*turb_scale
                     endif
#endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5d0)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5d0)*delta(2)+domnlo(2)

                  call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(1,2,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,1) = u + uflct(1,j,k)*turbSclX(x,y,z)*turb_scale
                     endif
#endif
                  endif

                  if (bc(1,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,2) = v + vflct(1,j,k)*turbSclY(x,y,z)*turb_scale
                     endif
#endif
                  endif

                  if (bc(1,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,3) = w + wflct(1,j,k)*turbSclZ(x,y,z)*turb_scale
                     endif
#endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5d0)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)

                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(2,1,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,1) = u + uflct(i,1,k)*turbSclX(x,y,z)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,2) = v + vflct(i,1,k)*turbSclY(x,y,z)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,1,3).eq.EXT_DIR) then
                        vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,3) = w + wflct(i,1,k)*turbSclZ(x,y,z)*turb_scale
                     endif

#endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5d0)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)

                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(2,2,1).eq.EXT_DIR) then
                        vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,1) = u + uflct(i,1,k)*turbSclX(x,y,z)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,2) = v + vflct(i,1,k)*turbSclY(x,y,z)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,3) = w + wflct(i,1,k)*turbSclZ(x,y,z)*turb_scale
                     endif
#endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5d0)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5d0)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
                  r = SQRT(x*x + y*y)
                  eta = EXP(-2*(Rfu-r)/.015d0)

                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(3,1,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 3) then
                        if(probtype.eq.4)then
                           vel(i,j,k,1) = u + uflct(i,j,1)*turb_scale
                        elseif(probtype.eq.5)then
                           vel(i,j,k,1) = u
                           if (r.le.Rfu) then
c                             vel(i,j,k,1) = u + uflct(i,j,1)*(2.7 + eta*4.25)
                              vel(i,j,k,1) = u + uflct(i,j,1)*(1.0d0 + eta*1.574d0)*factor
                              vel(i,j,k,1) = vel(i,j,k,1) * frcfct
                           endif
                        else
#ifdef SWIRL
                           vel(i,j,k,1) = u + uflct(i,j,1)
#else
                           r = sqrt(x**2+y**2)
                           cs = x/r
                           sg = y/r
                           scaler = turbSclR(r)
                           scalet = turbSclT(r)
                           vel(i,j,k,1) = u + turb_scale*turbSclX(x,y,z)*
     &                          ((scaler*cs**2+scalet*sg**2)*uflct(i,j,1) +
     &                          (scaler*sg*cs-scalet*sg*cs)*vflct(i,j,1))
#endif
                        endif
                     endif
#endif
                  endif

                  if (bc(3,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 3) then
                        if(probtype.eq.4)then
                           vel(i,j,k,2) = v + vflct(i,j,1)*turb_scale
                        elseif(probtype.eq.5)then
                           vel(i,j,k,2) = v
                           if (r.le.Rfu) then
c                             vel(i,j,k,2) = v + vflct(i,j,1)*(2.7 + eta*4.25)
                              vel(i,j,k,2) = v + vflct(i,j,1)*(1.0d0 + eta*1.574d0)*factor
                              vel(i,j,k,2) = vel(i,j,k,2) * frcfct
                           endif
                        else
#ifdef SWIRL
                           vel(i,j,k,2) = v + vflct(i,j,1)
#else
                           r = sqrt(x**2+y**2)
                           cs = x/r
                           sg = y/r
                           scaler = turbSclR(r)
                           scalet = turbSclT(r)
                           vel(i,j,k,2) = v + turb_scale*turbSclY(x,y,z)*
     &                          ((scaler*sg**2+scalet*cs**2)*vflct(i,j,1) +
     &                          (scaler*sg*cs- scalet*sg*cs)*uflct(i,j,1))
#endif
                        endif
                     endif
#endif
                  endif
                  
                  if (bc(3,1,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 3) then
                        if(probtype.eq.4)then
                           vel(i,j,k,3) = w + wflct(i,j,1)*turb_scale
                        elseif(probtype.eq.5)then
                           vel(i,j,k,3) = w
                           if (r.le.Rfu) then
c                             vel(i,j,k,3) = w + wflct(i,j,1)*(2.7 + eta*4.25)
                              vel(i,j,k,3) = w + wflct(i,j,1)*(1.0d0 + eta*1.574d0)*factor
                              vel(i,j,k,3) = vel(i,j,k,3) * frcfct
                           endif
                        else
#ifdef SWIRL
                           vel(i,j,k,3) = w + wflct(i,j,1)
#else 
                           vel(i,j,k,3) = w + wflct(i,j,1)*turbSclZ(x,y,z)*turb_scale
     &                          * anisotsc
#endif
                        endif
                     endif
#endif
                  endif
               enddo
            enddo
         enddo
      endif
      
      if (hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5d0)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5d0)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)

                  call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(3,2,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,1) = u + uflct(i,j,1)*turbSclX(x,y,z)*turb_scale
                     endif
#endif
                  endif

                  if (bc(3,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,2) = v + vflct(i,j,1)*turbSclY(x,y,z)*turb_scale
                     endif
#endif
                  endif

                  if (bc(3,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,3) = w + wflct(i,j,1)*turbSclZ(x,y,z)*turb_scale
                     endif
#endif
                  endif
               enddo
            enddo
         enddo
      endif

#if defined(BL_DO_FLCT)
      if (forceInflow) then
         deallocate(uflct)
         deallocate(vflct)
         deallocate(wflct)
      endif
#endif

      end

      subroutine FORT_XVELFILL (xvel,DIMS(xvel),domlo,domhi,delta,
     &                          xlo,time,bc)
      implicit none
      integer DIMDEC(xvel), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  xvel(DIMV(xvel))
      call bl_abort('should not be in XVELFILL')
      end

      subroutine FORT_YVELFILL (yvel,DIMS(yvel),domlo,domhi,delta,
     &                          xlo,time,bc)
      implicit none
      integer DIMDEC(yvel), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  yvel(DIMV(yvel))
      call bl_abort('should not be in YVELFILL')
      end

      subroutine FORT_ZVELFILL (zvel,DIMS(zvel),domlo,domhi,delta,
     &                          xlo,time,bc)
      implicit none
      integer DIMDEC(zvel), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  zvel(DIMV(zvel))
      call bl_abort('should not be in ZVELFILL')
      end

      subroutine FORT_ALLCHEMFILL (rhoY,DIMS(rhoY),domlo,domhi,delta,
     &                             xlo,time,bc)

      implicit none

#include "cdwrk.H"
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
            x = (float(i)+.5d0)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5d0)*delta(2)+domnlo(2)
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
            x = (float(i)+.5d0)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5d0)*delta(2)+domnlo(2)
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
            y = (float(j)+.5d0)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
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
            y = (float(j)+.5d0)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
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
            z = (float(k)+.5d0)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5d0)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
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
            z = (float(k)+.5d0)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5d0)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
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

#include "cdwrk.H"
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
            x = (float(i)+.5d0)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5d0)*delta(2)+domnlo(2)
                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5d0)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5d0)*delta(2)+domnlo(2)
                  call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5d0)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5d0)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5d0)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            z = (float(k)+.5d0)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5d0)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif    
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            z = (float(k)+.5d0)*delta(3)+domnlo(3)
            do j = lo(2),hi(2)
               y = (float(j)+.5d0)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
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

      end

      subroutine FORT_CALCDIVU(divu,DIMS(divu),plo,phi,xlo,dx,time)
c
c ... This routine is for user defined divu, for debugging purposes only.
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
               x = xlo(1) + (i-plo(1)+.5d0)*dx(1)
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
               divu(i,j,k) = time/.0125d0
            enddo
         enddo
      enddo
#endif
      do k=dlo(3),dhi(3)
         do j=dlo(2),dhi(2)
            do i = dlo(1), dhi(1)
               x = xlo(1) + (i-plo(1)+.5d0)*dx(1)
               divu(i,j,k) = merge(one,zero,x.GE.0.5d0)
            enddo
         enddo
      enddo

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
      REAL_T xx, yy, zz, r, dr, dz, blend
      
c     blend goes from 1 to 0 in r at wire_r, over a width of wire_w  
      blend(r) = one - half*(one+TANH((r-wire_r)/(fourth*wire_w)))

      do k = lo(3),hi(3)
         zz = float(k)*dx(3)+domnlo(3)
         do j = lo(2),hi(2)
            yy = float(j)*dx(2)+domnlo(2)
            do i = lo(1),hi(1)
               xx = float(i)*dx(1)+domnlo(1)

#if 0
c
c  Compute the distance to a disk at z=hotspot_z, radius hotspot_r and thickness
c     of hotspot_dz
c
               dr = SQRT(xx*xx+yy*yy) - hotspot_r
               if (dr .LT. zero) then
                  if (zz.GT.hotspot_z+hotspot_dz) then
                     r = zz - (hotspot_z+hotspot_dz)
                  else if (zz.LT.hotspot_z) then
                     r = hotspot_z - zz
                  else
                     r = zero
                  end if
               else if ((zz.GT.hotspot_z).AND.(zz.LT.hotspot_z+hotspot_dz)) then
                  r = MAX(dr,zero)
               else
                  dz = MIN(ABS(zz-hotspot_z-hotspot_dz),ABS(hotspot_z-zz))
                  r = SQRT(dr*dr+dz*dz)                  
               end if
               rad(i,j,k) = half*(one-tanh(two*(r-1.5d0*hotspot_w)/hotspot_w))
     &                                       *hotspot_htc*(hotspot_T - T(i,j,k))

#else
               r = SQRT((xx-wire_x)**2 + (zz-wire_z)**2 + MAX(zero,ABS(yy)-wire_y)**2)
               rad(i,j,k) = -blend(r)*wire_h*MAX(zero,(wire_T - T(i,j,k)))
#endif
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
c      
c     ::::: lets punt if not in domain interior
c
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

