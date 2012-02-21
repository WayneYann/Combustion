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

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

#define SETRHO_CALCT
#define SDIM 2

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
      REAL_T alpha,xsmb,vpmax,exp1 
      REAL_T  V_x, V_xx,beta,tauloc,hoft
      integer i1
      logical found_it

      do_control = 1

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
                  V_in = r2
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

      if (coft_old .lt. zero) coft_old = coft

      zbase_control = zbase_control + V_in*dt + dV_control*dt**2
      V_in_old = V_in
      V_in = V_in + dt*dV_control
      corr = 0.5d0

      slocal = half*(V_in_old + V_in) - (coft - coft_old)/(dt*scale_control)
      sest = (one - corr)*sest + corr*slocal
#if 0
      vslope = two*((cfix-coft)/(scale_control*tau_control) + sest - V_in)/tau_control

      V_new = V_in + dt*vslope
      dVmax = changeMax_control * one
#else



c     vpmax = max(V_in,V_in_old,one)*changeMax_control/dt
c     xsmb = -(cfix-coft)/scale_control

c     alpha = sqrt(abs(xsmb)/vpmax)

c     exp1 = dexp(-2.d0*tau_control/alpha)

c     V_new = sest+ (xsmb)/alpha*
c    &   (exp1*dexp(dt/alpha)-dexp(-dt/alpha))/
c    &      (1.d0+exp1)

      tauloc = MAX(tau_control,10.d0*dt)
      beta = (cfix-coft)/scale_control + sest * tauloc 
      V_x = 6.d0*(beta - tauloc*V_in - tauloc*(sest-V_in)/3.d0)/tauloc**2
      V_xx = 6.d0*(tauloc*(sest-V_in)/2.d0-(beta - tauloc*V_in))/tauloc**3

      V_new = V_in + dt * V_x + dt**2 * V_xx

      dVmax = changeMax_control * max(one,V_in)
#endif
c     V_new = MIN(MAX(V_new,V_in-dVmax),V_in+dVmax)
c     V_new = MIN(V_new,V_in+dVmax)
      V_new = MIN(MAX(V_new,V_in-2.d0*dVmax),V_in+changeMax_control * min(one,V_in))
      V_new = MAX(zero,V_new)

      tbase_control = time
      dV_control = (V_new - V_in)/dt

      hoft = coft / scale_control

      if (myproc.eq.0) then
         print *
         print *,'****************** control:'
         print *,'time,dt,coft,V_new:',time,dt,coft,V_new
         print *,'changeMax_control: ', changeMax_control
#if 1
         print *,'vslope,sest,cfix,slocal,V_x,V_xx:',vslope,sest,cfix,slocal,V_x,V_xx
#else
         print *,'alpha,sest,cfix,slocal:',alpha,sest,cfix,slocal
#endif
         print *,'hoft,coft_old,V_in_old:',hoft,coft_old,V_in_old
      endif

      coft_old = coft

      if (myproc.eq.0) then
         open(13,file=ac_hist_file,form='formatted',OPEN_APPEND)
         write(13,1000) step,time,V_in,tbase_control,zbase_control,
     &        dV_control,sest,coft_old
         close(13)
      endif
 1000 format(i7,8g26.18)
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

      integer i, j, n
      REAL_T Patm

     
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               scal(i,j,Trac) = zero
            end do
         end do
    

      Patm = pamb / 101325.0d0

      call FORT_RHOfromPTY(lo,hi,
     &     scal(ARG_L1(state),ARG_L2(state),Density), DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state),
     &     Patm)
      call FORT_HMIXfromTY(lo,hi,
     &     scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state))

     
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do n = 0,Nspec-1
                  scal(i,j,FirstSpec+n) = scal(i,j,FirstSpec+n)*scal(i,j,Density)
               enddo
               scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
            enddo
         enddo
    

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

#define RANDOM rand()
      REAL_T rand

#include "cdwrk.H"
#include "probdata.H"
#include "htdata.H"
#include "bc.H"
#include "visc.H"
#include "conp.H"
#if defined(BL_DO_FLCT)
#include "infl_frc.H"
#include "INFL_FORCE_F.H"
#endif

      REAL_T species(maxspec)
      integer lo(2), hi(2), n
      integer ifreq,zone,i,iH2,getZone
      character*(maxspnml) spName
      REAL_T x, y, time, u, v, rho, Yl(maxspec), Xl(maxspec), T, h, dx
      integer orient

      namelist /fortin/ vorterr, temperr, adverr, tempgrad, flametracval, probtype,
     &			max_temp_lev, max_vort_lev, splitx, splity, traceSpecVal,
     &                  xfrontw, yfrontw, refine_nozzle_x,
     &                  refine_nozzle_y, blobx, bloby,
     &                  blobr, xcen, v_strength,v_width,v_xcen,v_ycen,v_cl_x,
     &                  y_offset, v_blob_r, v_blob_T, v_blob_airfrac, V_in, V_co,
     &                  stTh, Rf, R_hot, R_hotBL, Twall, turb_scale, stBL,
     &                  max_trac_lev, max_nozzle_lev, pipeTh, pipeBL,
     &                  tV_in_l,tV_in_r,V_in_new,tV_co_l,tV_co_r,V_co_new,
     &                  T_stick, T_in,phi_in, thickFacTR, thickFacCH,nchemdiag,
     &                  pertmag
      namelist /heattransin/ pamb, dpdt_factor
#if defined(BL_DO_FLCT)
      namelist /flctin/ forceInflow, numInflPlanesStore, forceLo, forceHi,
     &     strmwse_dir, nCompInflow, flct_file
#endif
      namelist /control/ tau_control, sest, cfix, changeMax_control, h_control,
     &     zbase_control, tbase_control, pseudo_gravity
c
c      Build `probin' filename -- the name of file containing fortin namelist.
c
      integer maxlen, nCompFile, isioproc, ierr
      parameter (maxlen=256)

      character probin*(maxlen)

      call bl_pd_is_ioproc(isioproc)

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
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
      max_trac_lev = 100
      splitx = zero
      splity = zero
      traceSpecVal = 1.d-14
      xfrontw = fourth*half*(problo(1)+probhi(1))
      yfrontw = xfrontw
      refine_nozzle_y = zero
      refine_nozzle_x = zero
      max_nozzle_lev = 0
      blobx = -one
      bloby = -one
      blobr = -one
      xcen = problo(1)
      v_strength = zero
      v_width = -one
      v_xcen = zero
      v_ycen = zero
      v_cl_x = probhi(1)
      y_offset = zero
      pamb = 101325.d0
      dpdt_factor = 0.3d0
      v_blob_r = zero
      v_blob_T = 2000.d0
      v_blob_airfrac = zero
      V_in = zero
      V_co = zero
      stTh = zero
      Rf = 20*probhi(1)
      R_hot = half*probhi(1)
      R_hotBL = .1*probhi(1)
      Twall = 298.d0
      stBL = R_hotBL
      pipeTh = zero
      tV_in_l = 1.e6
      tV_in_r = 0.e6
      V_in_new = -100.d0
      tV_co_l = 1.e6
      tV_co_r = 0.e6
      V_co_new = -100.d0
      T_stick = -300.d0
      T_in = 300.d0
      phi_in = 0.d0
      zbase_control = 0.d0
      pertmag = 0.d0
      pseudo_gravity = 0
#if defined(BL_DO_FLCT)
c
c     Don't need to default 'nCompInflow' as it is block data'd to /3/
c
      forceInflow = .FALSE.
      numInflPlanesStore = -1
      forceLo = .TRUE.
      forceHi = .FALSE.
      strmwse_dir = FLCT_YVEL
      flct_file = ""
      turb_scale = 1
      nCompInFlow = 2
#endif

c     Initialize control variables
      do_control = 0
      tau_control = one
      sest = zero
      corr = one
      changeMax_control = .05
      coft_old = -one
      cfix = zero
      ac_hist_file = 'AC_History'
      h_control = -one
      nchemdiag = 1
      tbase_control = zero
      
      read(untin,fortin)

c     Initialize control variables that depend on fortin variables
      V_in_old = V_in
      
c      if (max_vort_lev.lt.0) max_vort_lev=max_temp_lev
      
      read(untin,heattransin)

#if defined(BL_DO_FLCT)
      read(untin,flctin)
#endif

      read(untin,control)
#ifndef GENGETFORCE
      if (pseudo_gravity.eq.1) then
         call bl_abort('Pseudo gravity requires -DGENGETFORCE')
      endif
#endif

   10 close(unit=untin)

#if defined(BL_DO_FLCT)
      if (forceInflow .eqv. .FALSE.) then
         forceLo = .FALSE.
         forceHi = .FALSE.
      else
         if (flct_file.ne."") then
#define FF_UNIT 20
            ierr = 0
            write(6,*) '...initializing turbulence, reading header info'
            open(FF_UNIT, file=trim(flct_file)//'/HDR',form='formatted',status='old',iostat=ierr)
            if (ierr .ne. 0) then
               call bl_abort('Problem opening file: ' // trim(flct_file) // '/HDR')
            end if
            call RD_SCL_FLCTHD(FF_UNIT,nCompFile,dimFile,probSizeFile,dxFile)
            close(FF_UNIT)
         endif
      endif
c      convVel = V_in
      convVel = one
#endif

c     Load domain dimensions into common, and set up boundary functions
      domnlo(1) = problo(1)
      domnlo(2) = problo(2)
      domnhi(1) = probhi(1)
      domnhi(2) = probhi(2)

      call setupbc()
      bcinit = .true.

      orient = 1
      x = half*(domnlo(1)+domnhi(1))
      y = domnlo(2)
      time = 0.d0
      call bcfunction(orient, x, y, time, u, v, rho, Yl, T, h, dx)
      call CKYTX(Yl,IWRK(ckbr),RWRK(ckbr),Xl)
      scale_control = Yl(fuelID)*rho*(domnhi(1)-domnlo(1))
      if (h_control .gt. zero) then
         cfix = scale_control * h_control
      else
         h_control = cfix / scale_control
      endif


      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
#if defined(BL_DO_FLCT)
         write(6,flctin)
#endif
         write(6,control)

         print *,'CONTROL STUFF:'
         print *,'  Fuel mass fraction:',Yl(fuelID)
         print *,'  Fuel mole fraction:',Xl(fuelID)
         print *,'  density at inlet:',rho
         print *,'  h_control:',h_control
         print *,'  scale_control:',scale_control
      end if

      end

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

#define M_FUEL  1
#define M_OXID  2
#define M_WALL  3
      
      subroutine setupbc()
      implicit none
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "htdata.H"
      
      REAL_T Patm, y, sum, yEval, pmfVals(maxspec+3), a, phi
      REAL_T X_bc(maxspec)
      REAL_T Yfu(maxspec), Yox(maxspec)
      REAL_T Xfu(maxspec), Xox(maxspec)
      integer i, zone, n, lo(SDIM), hi(SDIM), iN2, iH2, iO2, getZone
      integer isioproc
      character*(maxspnml) name
      data lo / 1, 1 /
      data hi / 1, 1 /
      
      call bl_pd_is_ioproc(isioproc)
      Patm = pamb / 101325.d0

      iN2 = -1
      do n = 1,Nspec
         call get_spec_name(name,n)
         if (name .eq. 'N2' ) iN2 = n
      enddo

c
c     For odd-numbered probdata, set M_FUEL mixture from the pmf file
c     For even-numbered probdata, set M_FUEL mixture using phi_in
c
      do zone=1,3

         do n=1,Nspec
            X_bc(n) = zero
         enddo        
         if (zone.eq.M_FUEL) then

            if ( MOD(probtype,2) .eq. 1) then

               yEval = (domnlo(2) - y_offset)*100.d0
               call pmf(yEval,yEval,pmfVals,n)      
               if ((isioproc.eq.1) .and. (n.ne.Nspec+3)) then
                  print *,'n=',n
                  print *,'Nspec=',Nspec
                  call bl_abort('SETUPBC: n .ne. Nspec+3')
               endif
               do n=1,Nspec
                  X_bc(n) = pmfVals(3+n)
               enddo        
               T_bc(zone) = pmfVals(1)
               u_bc(zone) = 0.d0
               v_bc(zone) = pmfVals(2)*1.d-2

            else

               phi = MIN( 1.d0, MAX( 0.d0, phi_in ) )
               call get_spec_name(name,fuelID)
               if (name .eq. 'H2' ) then
                  a = half
               else if (name .eq. 'CH4' ) then
                  a = two
               else if (name .eq. 'C3H8' ) then
                  a = five
               else
                  call bl_abort('Unknown fuel species')
               endif
               do n=1,Nspec
                  X_bc(n) = zero
               enddo        

               X_bc(oxidID) = one/(one + phi/a  + 0.79d0/0.21d0)
               X_bc(fuelID) = phi * X_bc(oxidID) / a
               X_bc(iN2) = 1.d0 - X_bc(oxidID) - X_bc(fuelID)
               T_bc(zone) = T_in
               u_bc(zone) = 0.d0
               v_bc(zone) = V_in

            endif

         else if (zone.eq.M_OXID) then
            
            X_bc(oxidID) = 0.21d0
            X_bc(iN2) = 1.d0 - X_bc(oxidID)
            T_bc(zone) = T_in
            u_bc(zone) = 0.d0
            v_bc(zone) = V_in
            
         else if (zone.eq.M_WALL) then
            
            X_bc(oxidID) = 0.21d0/0.79d0
            X_bc(iN2) = 1.d0 - Xfu(oxidID)
            T_bc(zone) = Twall
            u_bc(zone) = 0.d0
            v_bc(zone) = 0.d0
            
         else
            call bl_abort('Bad zone')
         endif

         call FORT_MOLETOMASS(lo,hi,
     &                        X_bc,        ARLIM(lo),ARLIM(hi),
     &                        Y_bc(0,zone),ARLIM(lo),ARLIM(hi))
         call FORT_RHOfromPTY(lo,hi,
     &                        rho_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        T_bc(zone),  ARLIM(lo),ARLIM(hi),
     &                        Y_bc(0,zone),ARLIM(lo),ARLIM(hi), Patm)
         call FORT_HMIXfromTY(lo,hi,
     &                        h_bc(zone),  ARLIM(lo),ARLIM(hi),
     &                        T_bc(zone),  ARLIM(lo),ARLIM(hi),
     &                        Y_bc(0,zone),ARLIM(lo),ARLIM(hi))
      enddo

c     Assume stick is in fuel, and compute reasonable rho,h there
      call FORT_RHOfromPTY(lo,hi,
     &                     rho_stick,     ARLIM(lo),ARLIM(hi),
     &                     T_stick,       ARLIM(lo),ARLIM(hi),
     &                     Y_bc(0,M_FUEL),ARLIM(lo),ARLIM(hi), Patm)
      call FORT_HMIXfromTY(lo,hi,
     &                     h_stick,       ARLIM(lo),ARLIM(hi),
     &                     T_stick,       ARLIM(lo),ARLIM(hi),
     &                     Y_bc(0,M_FUEL),ARLIM(lo),ARLIM(hi))
      end

         
c ::: -----------------------------------------------------------
      
      integer function getZone(x, y)
      implicit none
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      REAL_T x, y, eps, absx
      data eps / 1.d-8 /

      if (probtype.lt.5) then

         getZone = M_FUEL

         if (ABS(x-xcen) .gt. Rf + pipeTh) then
            getZone = M_OXID
         endif

         if (probtype.eq.4) then
            if ((x.lt.domnlo(1)+eps).or.(x.gt.domnhi(1)-eps)) then
               getZone = M_WALL
            endif
         endif
      else
         call bl_abort('Undefined problem type')
      endif
      end

c ::: -----------------------------------------------------------

      subroutine VAvg(V,x,y,dx,Msub)
      implicit none
#include "probdata.H"
#include "cdwrk.H"
#include "bc.H"
      REAL_T V(SDIM)
      REAL_T x,y,dx(SDIM)
      integer Msub

      REAL_T x1, incx
      integer i, zone, getZone

      incx = dx(1) / Msub
      V(1) = zero
      V(2) = zero
      do i=1,Msub
         x1 = ABS(x - half*dx(1) + (i-half)*incx)
         zone = getZone(x1,y)
         V(1) = V(1) + u_bc(zone)
         V(2) = V(2) + v_bc(zone)
      enddo
      V(1) = V(1)/Msub
      V(2) = V(2)/Msub
      end

c ::: -----------------------------------------------------------
      
      subroutine bcfunction(orient, x, y, time, u, v, rho, Yl, T, h, dx)
      implicit none
      integer orient
      REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(SDIM)
#include "cdwrk.H"
#include "htdata.H"
#include "bc.H"
#include "probdata.H"

      REAL_T Va(SDIM)
      integer n, getZone, zone, Msub
      parameter (Msub=1)

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if

      if (orient .lt. 4) then
         zone = getZone(x,y)

         u = u_bc(zone)
         v = v_bc(zone)
         rho = rho_bc(zone)
         do n = 0,Nspec-1
            Yl(n) = Y_bc(n,zone)
         end do
         T = T_bc(zone)
         h = h_bc(zone)

         if (do_control.eq.1) then

            if (zone .eq. M_FUEL) then
               u = 0
               v = V_in + (time-tbase_control)*dV_control
            endif
            
         else

c     Optionally, set time-dependent bc data
            if (zone .eq. M_FUEL) then
               if (tV_in_l.lt.tV_in_r) then
                  if (time .le. tV_in_l) then
                     v_bc(M_FUEL) = V_in
                  else if (time .ge. tV_in_r) then
                     v_bc(M_FUEL) = V_in_new
                  else
                     v_bc(M_FUEL) = V_in+(time-tV_in_l)*(V_in_new-V_in)/(tV_in_r-tV_in_l)
                  endif
               endif
               
            else if (zone .eq. M_OXID) then
               
               if (do_control.eq.0) then
                  if (tV_co_l .lt. tV_co_r) then 
                     if (time .le. tV_co_l) then
                        v_bc(M_OXID) = V_co
                     else if (time .ge. tV_co_r) then
                        v_bc(M_OXID) = V_co_new
                     else
                        v_bc(M_OXID) = V_co+(time-tV_co_l)*(V_co_new-V_co)/(tV_co_r-tV_co_l)
                     endif
                  endif
               endif
               
            endif
               
            call VAvg(Va,x,y,dx,Msub)         
            u = Va(1)
            v = Va(2)

         endif

c     If there is a flame-holder stick (gas mix will be fuel, and is loaded above
         if ((ABS(x).le..5*stTh) .and. (T_stick.gt.zero) .and. (stTh.gt.zero) ) then
            T = T_stick
            rho = rho_stick
            h = h_stick
         endif
               

      else
          print*, 'No boundary condition for orientation = ', orient
          call bl_abort(' ')
      end if

c      print*, x, v
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

#include "cdwrk.H"
#include "htdata.H"
#include "bc.H"
#include "probdata.H"

      integer i, j, n
      REAL_T  x, y
      integer slo(SDIM), shi(SDIM)
      character*(maxspnml) name
      REAL_T Patm
      integer spec, getZone, zone
      
      integer zL, zR, zH
      REAL_T eta, xi, xblend, mdotb, mdotL, mdotR, r1sq, r2sq, ylo, yhi
      REAL_T pmf_vals(maxspec+3), sum, Y_in(maxspec), rblob, theta
      REAL_T ypert, Lx, u, v, rho, T, h
      integer iN2, iO2
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1/
      data hi_chem /1,1/

      if ((Temp.gt.0).neqv.(RhoH.gt.0)) then
         call bl_abort('Need both Temp and RhoH, or neither')
      end if
      
      if ((Temp .LT. 0) .OR. (RhoH .LT. 0)) then
         call bl_abort('No IC''s for system without T, RhoH')
      endif

      iN2 = -1
      iO2 = -1
      do n = 1,Nspec
         call get_spec_name(name,n)
         if (name .eq. 'N2' ) iN2 = n
         if (name .eq. 'O2' ) iO2 = n
      end do

c     Set initial flame angle for probtype=2
      theta = Pi / 6.d0

c     Set vorticity into pressure field
c
      if (probtype.lt.4) then
         do j = ARG_L2(press),ARG_H2(press)
            y = float(j)*delta(2)+domnlo(2)
            do i = ARG_L1(press),ARG_H1(press)
               x = float(i)*delta(1)+domnlo(1)
               r1sq = (x-v_xcen)**2              + (y-v_ycen)**2
               r2sq = (x-(two*v_cl_x-v_xcen))**2 + (y-v_ycen)**2
               press(i,j) = v_strength*
     &              (dexp(-r1sq/v_width**2)-dexp(-r2sq/v_width**2)) 
            end do
         end do
         do j = lo(2), hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)         
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               
               ypert = 0.d0
               if (pertmag .gt. 0.d0) then
                  Lx = domnhi(1) - domnlo(1)
                  ypert = pertmag*(1.000 * sin(2*Pi*4*x/Lx)
     &                 + 1.023 * sin(2*Pi*2*(x-.004598)/Lx)
     &                 + 0.945 * sin(2*Pi*3*(x-.00712435)/Lx)
     &                 + 1.017 * sin(2*Pi*5*(x-.0033)/Lx)
     &                 + .982 * sin(2*Pi*5*(x-.014234)/Lx) )
               endif
               

               ylo = (y-half*delta(2)-y_offset+ypert)*100.d0
               yhi = (y+half*delta(2)-y_offset+ypert)*100.d0

               if (probtype.eq.2) then
                  ylo = ylo - MAX(zero,ABS(x)-half*stTh)*TAN(theta)
                  yhi = yhi - MAX(zero,ABS(x)-half*stTh)*TAN(theta)
               endif

               call pmf(ylo,yhi,pmf_vals,n)

               if (n.ne.Nspec+3) then
                  print *,'PMF file thinks n=',n
                  print *,'I think Nspec=',Nspec
                  call bl_abort('INITDATA: n .ne. Nspec+3, n:')
               endif

               call FORT_MOLETOMASS(lo_chem,hi_chem,
     &                              pmf_vals(4), ARLIM(lo_chem),ARLIM(hi_chem),
     &                               Y_in,       ARLIM(lo_chem),ARLIM(hi_chem))
               
               rblob = sqrt( (v_xcen-x)**2 + (v_ycen-y)**2 )
               
               scal(i,j,Trac) = fourth*(press(i,j  ) + press(i+1,j)
     &              +                   press(i,j+1) + press(i+1,j+1))

               if (v_blob_r .le. zero) then
                  scal(i,j,Temp) = pmf_vals(1)
                  do n = 1,Nspec
                     scal(i,j,firstSpec+n-1) = Y_in(n)
                  end do
                  if (iN2.gt.0) then
                     sum = zero
                     do n = 1,Nspec
                        if (n.ne.iN2) sum = sum+scal(i,j,firstSpec+n-1)
                     end do
                     scal(i,j,firstSpec+iN2-1) = one - sum
                  endif
               else
                  
c     Blend with air in blob, take air mix from bc(M_OXID)
                  eta = 0.5*(1.d0 + TANH((rblob - v_blob_r)/xfrontw))
                  
                  scal(i,j,Temp) = eta*pmf_vals(1) + (1.d0-eta)*v_blob_T
                  do n = 1,Nspec
                     scal(i,j,firstSpec+n-1) = eta*Y_in(n)
     &                    +(1.d0-eta)*(v_blob_airfrac*Y_bc(n-1,M_OXID)
     &                    +            (1.d0-v_blob_airfrac)*Y_in(n))
                  end do
                  
                  if (iN2.gt.0) then
                     sum = zero
                     do n = 1,Nspec
                        if (n.ne.iN2) sum = sum+scal(i,j,firstSpec+n-1)
                     end do
                     scal(i,j,firstSpec+iN2-1) = one - sum
                  endif
               end if
               
c     Cut off initial pmf at R_hot so that air along wall
               if (probtype.eq.2) then
                  eta = half*(one + TANH(two*(ABS(x)-R_hot)/R_hotBL))                  
                  scal(i,j,Temp) = (one - eta)*scal(i,j,Temp) + eta*T_in
                  do n = 1,Nspec
                     scal(i,j,firstSpec+n-1) =
     &                    (one - eta)*scal(i,j,firstSpec+n-1) 
     &                    +   eta * Y_bc(n-1,M_OXID)
                  end do 
               endif
               
               vel(i,j,1) = zero
               vel(i,j,2) = pmf_vals(2) / 1.d2
               
               if (probtype.eq.2) then
                  vel(i,j,2) = zero
               endif
            enddo
         enddo
      else
         do j = lo(2), hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)         
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               
#if 1
               eta = 0.5*(1.d0 + TANH((x - xcen)/xfrontw))
               do n = 0,Nspec-1
                  scal(i,j,firstSpec+n) =
     &                 (one - eta)*Y_bc(n,M_FUEL)
     &                 +   eta * Y_bc(n,M_OXID)
               end do 
               scal(i,j,Temp) = (1.d0-eta)*T_bc(M_FUEL) + eta*T_bc(M_OXID)
               scal(i,j,Trac) = 0.d0
               vel(i,j,1) = (1.d0-eta)*u_bc(M_FUEL) + eta*u_bc(M_OXID)
               vel(i,j,2) = (1.d0-eta)*v_bc(M_FUEL) + eta*v_bc(M_OXID)
#else
               y = domnlo(2)
               zone = getZone(x,y)
               call bcfunction(YLO, x, y, time, u, v, rho, Y_in, T, h, delta)
               do n = 0,Nspec-1
                  scal(i,j,firstSpec+n) = Y_in(n+1)
               end do 
               scal(i,j,Temp) = T
               scal(i,j,Trac) = 0.d0
               vel(i,j,1) = u
               vel(i,j,2) = v
#endif
            enddo
         enddo
      endif

      Patm = pamb / 101325.0
      call FORT_RHOfromPTY(lo,hi,
     &                     scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state),
     &                     scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),
     &                     scal(ARG_L1(state),ARG_L2(state),firstSpec),DIMS(state),
     &                     Patm)
      call FORT_HMIXfromTY(lo,hi,
     &                     scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state),
     &                     scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),
     &                     scal(ARG_L1(state),ARG_L2(state),firstSpec),DIMS(state))

c     Update typical values
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,Nspec-1
               typVal_Y(n+1) = MAX(typVal_Y(n+1),scal(i,j,firstSpec+n)) 
            enddo
            typVal_Density = MAX(scal(i,j,Density),typVal_Density)
            typVal_Temp    = MAX(scal(i,j,Temp),   typVal_Temp)
            typVal_Trac    = MAX(scal(i,j,Trac),   typVal_Trac)
            typVal_RhoH = MAX(ABS(scal(i,j,RhoH)*scal(i,j,Density)),typVal_RhoH)
            do n = 1,BL_SPACEDIM
               typVal_Vel  = MAX(ABS(vel(i,j,n)),typVal_Vel)
            enddo
         enddo
      enddo

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,Nspec-1
               scal(i,j,firstSpec+n) = scal(i,j,firstSpec+n)*scal(i,j,Density)
            enddo
            scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
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
c ::: lo,hi      => region of interest
c ::: domlo,hi   => index extent of problem domain
c ::: dx         => cell spacing
c ::: bc         => boundary condition flag (on orient)
c :::                   in BC_TYPES::physicalBndryTypes
c ::: problo     => phys loc of lower left corner of prob domain
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
      
#include "cdwrk.H"
#include "probdata.H"
#include "htdata.H"
      integer i, j, n, Tid, RHid, YSid, YEid, ys, ye
      integer getZone
      logical do_T, do_RH, do_Y
      REAL_T xl, xr, xh, y, ap

      if (probtype.eq.2) then
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
         if ((idir.EQ.1) .AND. (lo(2) .LE. domlo(2))
     &           .AND. (do_T .OR. do_RH .OR. do_Y) ) then
               
            j = lo(2)
            y = float(j)*dx(2)+domnlo(2)
            do i = lo(1), hi(1)
               
               xl = float(i)*dx(1)+domnlo(1) 
               xr = (float(i)+1.d0)*dx(1)+domnlo(1) 
               xh = 0.5d0*(xl+xr)
               
               if ( (getZone(xl,y).eq.M_FUEL) .OR.
     &              (getZone(xh,y).eq.M_FUEL) .OR.
     &              (getZone(xr,y).eq.M_FUEL) ) then
                  
#if 0
                  ap = zero
                  if (xl.ge.-half*stTh) then
                     ap = ap + (MIN(Rf,xr)-MAX(xl,half*stTh))/dx(1)
                  else if (xr.le.half*stTh) then
                     ap = ap + (MIN(xr,-half*stTh)-MAX(-Rf,xl))/dx(1)
                  else
                     call bl_abort("cell spans zero")
                  endif
#else
                  ap = zero
#endif

#if 0
                  if (do_T)  diff(i,j,Tid ) = diff(i,j,Tid )*ap
                  if (do_RH) diff(i,j,RHid) = diff(i,j,RHid)*ap
#else
c     Conducting stick: only zero conduction if not on stick
                  if ( (T_stick .le. zero)  .or. (
     &                 (ABS(xl).gt..5*stTh) .and.
     &                 (ABS(xh).gt..5*stTh) .and.
     &                 (ABS(xr).gt..5*stTh) ) ) then

                     if (do_T)  diff(i,j,Tid ) = diff(i,j,Tid )*ap
                     if (do_RH) diff(i,j,RHid) = diff(i,j,RHid)*ap

                  endif
#endif
                  if (do_Y) then
                     do n=ys,ye
                        diff(i,j,n) = diff(i,j,n)*ap
                     enddo
                  endif
               endif
            enddo
         endif
      endif
      
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

      integer   i, j
      REAL_T    x, y
      logical   in_refine_zone

#include "probdata.H"

      if (level .lt. max_trac_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j) = merge(set,tag(i,j),ftrac(i,j,1).gt.flametracval)
            enddo
         enddo
      endif

      if (level .lt. max_nozzle_lev) then
         do j = lo(2), hi(2)
            y = (float(j)+.5)*dx(2)+problo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*dx(1)+problo(1)
               in_refine_zone = (y - problo(2)) .le. refine_nozzle_y
               in_refine_zone = in_refine_zone .and.
     &              (ABS(x-xcen) .le. refine_nozzle_x)
               tag(i,j) = merge(set,tag(i,j),in_refine_zone)
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

      if (time .eq. zero) then
         call FORT_MVERROR(tag,DIMS(tag),set,clear,
     &                     adv,DIMS(adv),lo,hi,nvar,
     &                     domlo,domhi,delta,xlo,
     &			   problo,time,level)
      end if

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

      REAL_T    x, y, ax, ay, aerr
      integer   i, j, ng

#include "probdata.H"

      ng = min(ARG_H1(temp)-hi(1),ARG_H2(temp)-hi(2),lo(1)-ARG_L1(temp),lo(2)-ARG_L2(temp))

      if (ng .lt. 1) then
         print*, "TEMPERR cannot compute gradient, ng = ",ng
         call bl_abort(' ')
      endif
c
c     refine where there is temperature gradient
c
      if (level .lt. max_temp_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
c              ax = abs(temperature(i+1,j,1) - temperature(i-1,j,1))
c              ay = abs(temperature(i,j+1,1) - temperature(i,j-1,1))
c              aerr = max(ax,ay)
c              tag(i,j) = merge(set,tag(i,j),aerr.ge.tempgrad)
               tag(i,j) = merge(set,tag(i,j),temperature(i,j,1).le.tempgrad)
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

      REAL_T    x, y
      integer   i, j
      REAL_T    max_vort, tagval, v, r1sq, r2sq

#include "probdata.H"

      if (time.eq.zero .and. level.lt.max_vort_lev) then
         do j = lo(2), hi(2)
            y = float(j)*dx(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = float(i)*dx(1)+domnlo(1)
               r1sq = (x-v_xcen)**2              + (y-v_ycen)**2
               r2sq = (x-(two*v_cl_x-v_xcen))**2 + (y-v_ycen)**2
               v = v_strength*
     &              (dexp(-r1sq/v_width**2)-dexp(-r2sq/v_width**2)) 
               tag(i,j) = merge(set,tag(i,j),ABS(v).ge.vorterr)
            end do
         end do
      endif

      if (level .lt. max_vort_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j) = merge(set,tag(i,j),ABS(vort(i,j,1)).ge.vorterr*2.d0**level)
            end do
         end do
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
      
      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx, xhi(SDIM)
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(den)
      hi(1) = ARG_H1(den)
      lo(2) = ARG_L2(den)
      hi(2) = ARG_H2(den)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      
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
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
               den(i,j) = rho
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
               den(i,j) = rho
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
               den(i,j) = rho
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
               den(i,j) = rho
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

      integer    i, j
      integer    ilo, ihi, jlo, jhi
      REAL_T     y

#include "probdata.H"

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(adv)
      hi(1) = ARG_H1(adv)
      lo(2) = ARG_L2(adv)
      hi(2) = ARG_H2(adv)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))

      call filcc (adv,DIMS(adv),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
#if 0
        stop 'ADVFILL: should not be here for pipe problem'
#endif
           do i = lo(1), domlo(1)-1
             do j = lo(2), hi(2)
               adv(i,j) = zero
             enddo
           enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
#if 0
        stop 'ADVFILL: should not be here for pipe problem'
#endif
           do i = domhi(1)+1, hi(1)
             do j = lo(2), hi(2)
               adv(i,j) = zero
             enddo
           enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
c                               inflow for burner in a can (bic, biac)

           do j = lo(2), domlo(2)-1
              do i = lo(1), hi(1)
                 adv(i,j) = zero
              enddo
           enddo
      endif            

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
#if 0
        stop 'ADVFILL: should not be here for pipe problem'
#endif
           do j = domhi(2)+1, hi(2)
              do i = lo(1), hi(1)
                 adv(i,j) = zero
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
     &                          xlo,time,bc)
      implicit none
      integer DIMDEC(temp), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  temp(DIMV(temp))

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx, xhi(SDIM)
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(temp)
      hi(1) = ARG_H1(temp)
      lo(2) = ARG_L2(temp)
      hi(2) = ARG_H2(temp)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      
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
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
               temp(i,j) = T
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
               temp(i,j) = T
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
               temp(i,j) = T
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
               temp(i,j) = T
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
     &                          xlo,time,bc,id )
      implicit none
      integer DIMDEC(rhoh), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM), id
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoh(DIMV(rhoh))

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx, xhi(SDIM)
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(rhoh)
      hi(1) = ARG_H1(rhoh)
      lo(2) = ARG_L2(rhoh)
      hi(2) = ARG_H2(rhoh)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      
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
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif
      
      end
c
c Fill x & y velocity at once.
c
      subroutine FORT_VELFILL (vel,DIMS(vel),domlo,domhi,delta,
     &                         xlo,time,bc)

      implicit none
      integer DIMDEC(vel), bc(SDIM,2,SDIM)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  vel(DIMV(vel),SDIM)

      call FORT_XVELFILL (vel(ARG_L1(vel),ARG_L2(vel),1),
     & DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,1))

      call FORT_YVELFILL (vel(ARG_L1(vel),ARG_L2(vel),2),
     & DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,2))

      end

c
c Fill all chem species at once
c
      subroutine FORT_ALLCHEMFILL (rhoY,DIMS(rhoY),domlo,domhi,delta,
     &                             xlo,time,bc)

      implicit none
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"

      integer DIMDEC(rhoY), bc(SDIM,2,Nspec)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoY(DIMV(rhoY),Nspec)

      integer n
      
      do n=1,Nspec
         call FORT_CHEMFILL (rhoY(ARG_L1(rhoY),ARG_L2(rhoY),n),
     &        DIMS(rhoY),domlo,domhi,delta,xlo,time,bc(1,1,n),n-1)
      enddo
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

      REAL_T function turbSclX(x,y)
#include "probdata.H"
      REAL_T x, y, shapet, eta
      if ((probtype.eq.1).or.(probtype.eq.3)) then
         turbSclX = one
      else
         eta = ABS(x)/Rf
         turbSclX = shapet(eta)
         if (stBL.gt.zero) then
            turbSclX = turbSclX*TANH(4.d0*MAX(0.d0,ABS(x)-0.5d0*stTh)/stBL)
         end if
      end if
      end

      REAL_T function turbSclY(x,y)
#include "probdata.H"
      REAL_T x, y, shapen, eta
      turbSclY = zero
      if ((probtype.eq.1).or.(probtype.eq.3)) then
         turbSclY = one
      else
         eta = ABS(x)/Rf
         if(eta.lt.one)then
            turbSclY = 2.d0 -1.d0 * tanh((1-eta)/.0428222042d0)
         endif
         if (stBL.gt.zero) then
            turbSclY = turbSclY*TANH(4.d0*MAX(0.d0,ABS(x)-0.5d0*stTh)/stBL)
         end if
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

#if defined(BL_DO_FLCT)
#include "INFL_FORCE_F.H"
      REAL_T turbSclX, t_flct, dt_flct
      integer DIMDEC(uflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T uflct(:,:)
      allocatable uflct
#endif
      
      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx, xhi(SDIM)
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(xvel)
      hi(1) = ARG_H1(xvel)
      lo(2) = ARG_L2(xvel)
      hi(2) = ARG_H2(xvel)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      
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
         dt_flct = time - tbase_control
         t_flct = zbase_control + V_in*dt_flct + dV_control*dt_flct**2
c         write(*,*) 'time,dt_flct,t_flct,tbase_control:',time,dt_flct,t_flct,tbase_control
         call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo, delta, t_flct,
     $                  bc, domnlo, domnhi)
      endif
#endif

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
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
               xvel(i,j) = u
#if defined(BL_DO_FLCT)
               if (forceLo .and. strmwse_dir .eq. 1) then
                  xvel(i,j) = xvel(i,j) + uflct(1,j)*turb_scale
               endif
#endif
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
               xvel(i,j) = u
#if defined(BL_DO_FLCT)
               if (forceHi .and. strmwse_dir .eq. 1) then
                  xvel(i,j) = xvel(i,j) + uflct(1,j)*turb_scale
               endif
#endif
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
               xvel(i,j) = u
#if defined(BL_DO_FLCT)
               if (forceLo .and. strmwse_dir .eq. 2) then
                  xvel(i,j) = xvel(i,j) + uflct(i,1)*turb_scale*turbSclX(x,y)
               endif
#endif
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
               xvel(i,j) = u
#if defined(BL_DO_FLCT)
               if (forceHi .and. strmwse_dir .eq. 2) then
                  xvel(i,j) = xvel(i,j) + uflct(i,1)*turb_scale
               endif
#endif
            enddo
         enddo
      endif
      
#if defined(BL_DO_FLCT)
      if (forceInflow) then
         deallocate(uflct)
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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      
#if defined(BL_DO_FLCT)
#include "INFL_FORCE_F.H"
      REAL_T turbSclY, t_flct, dt_flct
      integer DIMDEC(vflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T vflct(:,:)
      allocatable vflct
#endif

      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx, xhi(SDIM)
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(yvel)
      hi(1) = ARG_H1(yvel)
      lo(2) = ARG_L2(yvel)
      hi(2) = ARG_H2(yvel)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      
#if defined(BL_DO_FLCT)
      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         enddo
         loFlctArray(strmwse_dir) = 1
         hiFlctArray(strmwse_dir) = 1
         call SET_ARGS(DIMS(vflct), loFlctArray, hiFlctArray)
         allocate(vflct(DIMV(vflct)))
         dt_flct = time - tbase_control
         t_flct = zbase_control + V_in*dt_flct + dV_control*dt_flct**2
         call INFL_FILL(FLCT_YVEL, DIMS(vflct), vflct, xlo, delta, t_flct,
     $                  bc, domnlo, domnhi)
      endif
#endif

      call filcc (yvel,DIMS(yvel),domlo,domhi,delta,xlo,bc)

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
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
               yvel(i,j) = v
#if defined(BL_DO_FLCT)
               if (forceLo .and. strmwse_dir .eq. 1) then
                  yvel(i,j) = yvel(i,j) + vflct(1,j)*turb_scale
               endif
#endif
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
               yvel(i,j) = v
#if defined(BL_DO_FLCT)
               if (forceHi .and. strmwse_dir .eq. 1) then
                  yvel(i,j) = yvel(i,j) + vflct(1,j)*turb_scale
               endif
#endif
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
               yvel(i,j) = v
#if defined(BL_DO_FLCT)
               if (forceLo .and. strmwse_dir .eq. 2) then
                  yvel(i,j) = yvel(i,j) + vflct(i,1)*turb_scale*turbSclY(x,y)
               endif
#endif
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
               yvel(i,j) = v
#if defined(BL_DO_FLCT)
               if (forceHi .and. strmwse_dir .eq. 2) then
                  yvel(i,j) = yvel(i,j) + vflct(i,1)*turb_scale
               endif
#endif
            enddo
         enddo
      endif
      
#if defined(BL_DO_FLCT)
      if (forceInflow) then
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
     &                          xlo,time,bc,id )
      implicit none
      integer DIMDEC(rhoY), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM), id
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoY(DIMV(rhoY))

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      
      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx, xhi(SDIM)
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(rhoY)
      hi(1) = ARG_H1(rhoY)
      lo(2) = ARG_L2(rhoY)
      hi(2) = ARG_H2(rhoY)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      
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
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
               rhoY(i,j) = rho*Yl(id)
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

      integer    i, j
      integer    ilo, ihi, jlo, jhi
      logical    fix_xlo, fix_xhi, fix_ylo, fix_yhi
      logical    per_xlo, per_xhi, per_ylo, per_yhi

      fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
      per_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .eq. INT_DIR)
      fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
      per_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .eq. INT_DIR)
      fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
      per_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .eq. INT_DIR)
      fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)
      per_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .eq. INT_DIR)

      ilo = max(ARG_L1(p),domlo(1))
      ihi = min(ARG_H1(p),domhi(1))
      jlo = max(ARG_L2(p),domlo(2))
      jhi = min(ARG_H2(p),domhi(2))
c
c     ::::: left side
c
      if (fix_xlo) then
         do i = ARG_L1(p), domlo(1)-1
            do j = jlo,jhi
               p(i,j) = p(ilo,j)
            end do
         end do
         if (fix_ylo) then
            do i = ARG_L1(p), domlo(1)-1
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ilo,jlo)
               end do
            end do
         else if (per_ylo) then
            do i = ARG_L1(p), domlo(1)-1
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ilo,j)
               end do
            end do
         end if
         if (fix_yhi) then
            do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ilo,jhi)
               end do
            end do
         else if (per_yhi) then
            do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ilo,j)
               end do
            end do
         end if
      end if
c
c     ::::: right side
c
      if (fix_xhi) then
         do i = domhi(1)+1, ARG_H1(p)
            do j = jlo,jhi
               p(i,j) = p(ihi,j)
            end do
	 end do
	 if (fix_ylo) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ihi,jlo)
               end do
	    end do
	 else if (per_ylo) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ihi,j)
               end do
	    end do
         end if
	 if (fix_yhi) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ihi,jhi)
               end do
	    end do
	 else if (per_yhi) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ihi,j)
               end do
	    end do
         end if
      end if
      
      if (fix_ylo) then
         do j = ARG_L2(p), domlo(2)-1
            do i = ilo, ihi
               p(i,j) = p(i,jlo)
            end do
	 end do
	 if (per_xlo) then
          do j = ARG_L2(p), domlo(2)-1
               do i = ARG_L1(p), domlo(1)-1
                  p(i,j) = p(i,jlo)
               end do
	    end do
         end if
	 if (per_xhi) then
           do j = ARG_L2(p), domlo(2)-1
               do i = domhi(1)+1, ARG_H1(p)
                  p(i,j) = p(i,jlo)
               end do
	    end do
         end if
      end if

      if (fix_yhi) then
         do j = domhi(2)+1, ARG_H2(p)
            do i = ilo, ihi
               p(i,j) = p(i,jhi)
            end do
	 end do
	 if (per_xlo) then
	    do j = domhi(2)+1, ARG_H2(p)
               do i = ARG_L1(p), domlo(1)-1
                  p(i,j) = p(i,jhi)
               end do
	    end do
         end if
	 if (per_xhi) then
	    do j = domhi(2)+1, ARG_H2(p)
               do i = domhi(1)+1, ARG_H1(p)
                  p(i,j) = p(i,jhi)
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

      integer    i,j
      REAL_T     x,time

      dlo(1) = divu_l1
      dlo(2) = divu_l2
      dhi(1) = divu_h1
      dhi(2) = divu_h2

#if 0
      do j=dlo(2),dhi(2)
        do i = dlo(1), dhi(1)
          x = xlo(1) + (i-plo(1)+.5)*dx(1)
          if (x.le..005) then
            divu(i,j) = zero
          else
c            divu(i,j) = one
c            divu(i,j) = min(one,time/.05)
c            divu(i,j) = min(two,one+time/.05)
            divu(i,j) = time
          endif
        enddo
      enddo
#endif
#if 0
      do j=dlo(2),dhi(2)
        do i = dlo(1), dhi(1)
            divu(i,j) = time/.0125
        enddo
      enddo
#endif
      do j=dlo(2),dhi(2)
        do i = dlo(1), dhi(1)
           x = xlo(1) + (i-plo(1)+.5)*dx(1)
           divu(i,j) = merge(one,zero,x.GE.0.5)
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

      integer i, j
      
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rad(i,j) = zero
         end do
      end do
      end

c
c
c ::: -----------------------------------------------------------
c
c     This routine add the forcing terms to the momentum equation
c
      subroutine FORT_MAKEFORCE(time,force,rho,
     &                          DIMS(istate),DIMS(state),
     &                          dx,xlo,xhi,gravity,scomp,ncomp)

      implicit none

      integer    DIMDEC(state)
      integer    DIMDEC(istate)
      integer    scomp, ncomp
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     force  (DIMV(istate),scomp+1:scomp+ncomp)
      REAL_T     rho    (DIMV(state))
      REAL_T     gravity

#include "probdata.H"
#include "cdwrk.H"
#include "bc.H"

c
c     ::::: local variables
c
      integer i, j, n
      integer ilo, jlo
      integer ihi, jhi
      integer a2, a3, a4, a5
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  sga, cga
      integer isioproc
      integer nXvel, nYvel, nRho, nTrac

      call bl_pd_is_ioproc(isioproc)

      if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
         write(*,*) "pseudo_gravity::dV_control = ",dV_control
      endif

      hx = dx(1)
      hy = dx(2)

      ilo = istate_l1
      jlo = istate_l2
      ihi = istate_h1
      jhi = istate_h2

c     Assumes components are in the following order
      nXvel = 1
      nYvel = 2
      nRho  = 3
      nTrac = 4

      if (scomp.eq.0) then
         if (abs(gravity).gt.0.0001) then
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nXvel) = zero
                  force(i,j,nYvel) = gravity*rho(i,j)
               enddo
            enddo
c     else to zero
         else
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nXvel) = zero
                  force(i,j,nYvel) = zero
               enddo
            enddo
         endif
c     Let's add the pseudo gravity afterwards...
         if (pseudo_gravity.eq.1) then
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nYvel) = force(i,j,nYvel) + dV_control*rho(i,j)
               enddo
            enddo
         endif
c     End of velocity forcing
      endif
      
      if ((scomp+ncomp).gt.BL_SPACEDIM) then
c     Scalar forcing
         do n = max(scomp+1,nRho), scomp+ncomp
            if (n.eq.nRho) then
c     Density
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else if (n.eq.nTrac) then
c     Tracer
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else
c     Other scalar
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            endif
         enddo
      endif

      end

