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

      if (coft_old.lt.zero) coft_old = coft

      if (dt.gt.zero) then
         zbase_control = zbase_control + V_in*dt + dV_control*dt**2
         V_in_old = V_in
         V_in = V_in + dt*dV_control

         slocal = half*(V_in_old + V_in) - (coft - coft_old)/(dt*scale_control)
         sest = (one - corr)*sest + corr*slocal
         vslope = two*((cfix-coft)/(scale_control*tau_control) + sest - V_in)/tau_control
         
         V_new = V_in + dt*vslope
         dVmax = changeMax_control *  MAX(V_in_old,one)
         V_new = MIN(MAX(V_new,V_in-dVmax),V_in+dVmax)
         V_new = MAX(zero,V_new)
         
         tbase_control = time
         dV_control = (V_new - V_in)/dt
         
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
            open(13,file=ac_hist_file,form='formatted',position='append')
            write(13,1000) step,time,V_in,tbase_control,zbase_control,
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
#include "cdwrk.H"
#include "htdata.H"
#include "bc.H"
#if defined(BL_DO_FLCT)
#include "INFL_FORCE_F.H"
#endif
#include "visc.H"
#include "conp.H"

#ifdef DO_LMC_FORCE
#include "forcedata.H"
#endif

      integer lo(SDIM), hi(SDIM), i, j, zone, iH2, getZone
      character*(maxspnml) spName
      REAL_T rand

      namelist /fortin/ vorterr, temperr, adverr, tempgrad, 
     &                  flametracval, probtype,
     &		        max_temp_lev, max_vort_lev, max_trac_lev,
     &                  traceSpecVal,phi_in,T_in,
     $                  turb_scale, V_in, 
     $                  zstandoff, pertmag, nchemdiag
      namelist /heattransin/ pamb, dpdt_factor
#if defined(BL_DO_FLCT)
      namelist /flctin/ tstart_turb, forceInflow, numInflPlanesStore, forceLo, forceHi,
     &     strmwse_dir, nCompInflow, flct_file
#endif
      namelist /control/ tau_control, sest, cfix, changeMax_control, h_control, 
     &     zbase_control, pseudo_gravity
#ifdef DO_LMC_FORCE
      REAL_T  twicePi, kxd, kyd, kzd
      REAL_T  thetaTmp, phiTmp
      REAL_T  cosThetaTmp, cosPhiTmp
      REAL_T  sinThetaTmp, sinPhiTmp
      REAL_T  px, py, pz, mp2, Ekh
      integer kx, ky, kz, mode_count, reduced_mode_count
      integer xstep, ystep, zstep

      REAL_T  Lx, Ly, Lz, Lmin, rn 
      REAL_T  kappa, kappaMax, freqMin, freqMax, freqDiff, pdk

      namelist /fortin/ nmodes, nxmodes, nymodes, nzmodes, mode_start, hack_lz,
     &                  forcing_type, spectrum_type, time_offset, forcing_twice_wavelength,
     &                  forcing_xlength, forcing_ylength, forcing_zlength,
     &                  forcing_time_scale_min, forcing_time_scale_max, 
     &                  force_scale, forcing_epsilon, blrandseed,
     &                  use_rho_in_forcing, do_mode_division,
     &                  div_free_force, moderate_zero_modes
#endif

c
c      Build `probin' filename -- the name of file containing fortin namelist.
c
      integer maxlen, nCompFile, isioproc, ierr
      parameter (maxlen=256)

      character probin*(maxlen)
      REAL_T y, yEval, pmf_vals(maxspec+3)
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
      max_trac_lev = 100
      traceSpecVal = 1.d-10
      pamb = 101325.d0
      dpdt_factor = 0.3d0

#if defined(BL_DO_FLCT)
c     add_turb = .FALSE.
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
      turb_scale = 1
#endif

      zbase_control = 0.d0

c     Note: for setup with no coflow, set Ro=Rf+wallth
      zstandoff = zero
      pertmag = 0.d0

c     Initialize control variables
      tau_control = one
      sest = zero
      corr = one
      changeMax_control = .05
      coft_old = -one
      cfix = zero
      ac_hist_file = 'AC_History'
      h_control = -one
      nchemdiag = 1
      dV_control = zero
      tbase_control = zero
      h_control = -one
      pseudo_gravity = 0
      nchemdiag = 1

      read(untin,fortin)
      
c     Initialize control variables that depend on fortin variables
      V_in_old = V_in

      if (max_vort_lev.lt.0) max_vort_lev=max_temp_lev
      
      read(untin,heattransin)
 
#if defined(BL_DO_FLCT)
      read(untin,flctin)
#endif
 
      read(untin,control)

   10 close(unit=untin)

#if defined(BL_DO_FLCT)
      if (forceInflow .eqv. .FALSE.) then
         forceLo = .FALSE.
         forceHi = .FALSE.
      else
         if (flct_file.ne.'') then
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
      convVel = V_in
      if (probtype.ge.3) then
         convVel = one
      endif
#endif
c     Load domain dimensions into common, and set up boundary functions
      domnlo(1) = problo(1)
      domnlo(2) = problo(2)
      domnlo(3) = problo(3)
      domnhi(1) = probhi(1)
      domnhi(2) = probhi(2)
      domnhi(3) = probhi(3)
      
      call setupbc()
      bcinit = .true.
      
      iH2 = -1
      do n = 1,Nspec
         call get_spec_name(spName,n)
         if (spName .eq. 'H2') iH2 = n-1
      end do
      zone = getZone(half*(domnlo(1)+domnhi(1)),half*(domnlo(2)+domnhi(2)),domnlo(3))
      scale_control = Y_bc(iH2,zone)*rho_bc(zone)*(domnhi(1)-domnlo(1))*(domnhi(2)-domnlo(2))
      if (h_control .gt. zero) then
         cfix = scale_control * h_control
      endif

#ifdef DO_LMC_FORCE
      if ( probtype.eq.3 ) then

         if (isioproc .eq. 1) then
            write (*,*) "Initialising random number generator..."
         endif
         
         twicePi = two*Pi
         
         if (blrandseed.gt.0) then
            call blutilinitrand(blrandseed)
            call blutilrand(rn)
            call blutilinitrand(blrandseed)
            if (isioproc .eq. 1) then
               write (*,*) "blrandseed = ",blrandseed
               write (*,*) "first random number = ",rn
            endif
         else
            call blutilinitrand(111397)
         endif

         Lx = domnhi(1)-domnlo(1)
         Ly = domnhi(2)-domnlo(2)
         Lz = domnhi(3)-domnlo(3)
         
         if (hack_lz.eq.1) then 
            Lz = Lz/two
         endif
         
         if (isioproc .eq. 1) then
            write(*,*) "Lx = ",Lx
            write(*,*) "Ly = ",Ly 
            write(*,*) "Lz = ",Lz
         endif
         
         Lmin = min(Lx,Ly,Lz)
         kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
         nxmodes = nmodes*int(0.5+Lx/Lmin)
         nymodes = nmodes*int(0.5+Ly/Lmin)
         nzmodes = nmodes*int(0.5+Lz/Lmin)
         if (isioproc .eq. 1) then
            write(*,*) "Lmin = ",Lmin
            write(*,*) "kappaMax = ",kappaMax
            write(*,*) "nxmodes = ",nxmodes
            write(*,*) "nymodes = ",nymodes
            write(*,*) "nzmodes = ",nzmodes
         endif
         
         if (forcing_time_scale_min.eq.zero) then
            forcing_time_scale_min = half
         endif
         if (forcing_time_scale_max.eq.zero) then
            forcing_time_scale_max = one
         endif

         freqMin = one/forcing_time_scale_max
         freqMax = one/forcing_time_scale_min
         freqDiff= freqMax-freqMin
         
         if (isioproc .eq. 1) then
            write(*,*) "forcing_time_scale_min = ",forcing_time_scale_min
            write(*,*) "forcing_time_scale_max = ",forcing_time_scale_max
            write(*,*) "freqMin = ",freqMin
            write(*,*) "freqMax = ",freqMax
            write(*,*) "freqDiff = ",freqDiff
         endif
         
         mode_count = 0

         xstep = int(Lx/Lmin+0.5)
         ystep = int(Ly/Lmin+0.5)
         zstep = int(Lz/Lmin+0.5)
         if (isioproc .eq. 1) then
            write (*,*) "Mode step ",xstep, ystep, zstep
         endif
         
         do kz = mode_start*zstep, nzmodes, zstep
            kzd = dfloat(kz)
            do ky = mode_start*ystep, nymodes, ystep
               kyd = dfloat(ky)
               do kx = mode_start*xstep, nxmodes, xstep
                  kxd = dfloat(kx)
                  
                  kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                  
                  if (kappa.le.kappaMax) then
                     call blutilrand(rn)
                     FTX(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                     call blutilrand(rn)
                     FTY(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                     call blutilrand(rn)
                     FTZ(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
c     Translation angles, theta=0..2Pi and phi=0..Pi
                     call blutilrand(rn)
                     TAT(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     TAP(kx,ky,kz) = rn*Pi
c     Phases
                     call blutilrand(rn)
                     FPX(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPY(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPZ(kx,ky,kz) = rn*twicePi
                     if (div_free_force.eq.1) then
                        call blutilrand(rn)
                        FPXX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPXY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPXZ(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYZ(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZZ(kx,ky,kz) = rn*twicePi
                     endif
c     Amplitudes (alpha)
                     call blutilrand(rn)
                     thetaTmp      = rn*twicePi
                     cosThetaTmp   = cos(thetaTmp)
                     sinThetaTmp   = sin(thetaTmp)
                     call blutilrand(rn)
                     phiTmp        = rn*Pi
                     cosPhiTmp     = cos(phiTmp)
                     sinPhiTmp     = sin(phiTmp)

                     px = cosThetaTmp * sinPhiTmp
                     py = sinThetaTmp * sinPhiTmp
                     pz =               cosPhiTmp

                     mp2           = px*px + py*py + pz*pz
                     if (kappa .lt. 0.000001) then
                        if (isioproc .eq. 1) then
                           write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                        endif
                        FAX(kx,ky,kz) = zero
                        FAY(kx,ky,kz) = zero
                        FAZ(kx,ky,kz) = zero
                     else
c     Count modes that contribute
                        mode_count = mode_count + 1
c     Set amplitudes
                        if (spectrum_type.eq.1) then
                           Ekh        = one / kappa
                        else if (spectrum_type.eq.2) then
                           Ekh        = one / (kappa*kappa)
                        else
                           Ekh        = one
                        endif
                        if (div_free_force.eq.1) then
                           Ekh = Ekh / kappa
                        endif
                        if (moderate_zero_modes.eq.1) then
                           if (kx.eq.0) Ekh = Ekh / two
                           if (ky.eq.0) Ekh = Ekh / two
                           if (kz.eq.0) Ekh = Ekh / two
                        endif
                        if (force_scale.gt.zero) then
                           FAX(kx,ky,kz) = force_scale * px * Ekh / mp2
                           FAY(kx,ky,kz) = force_scale * py * Ekh / mp2
                           FAZ(kx,ky,kz) = force_scale * pz * Ekh / mp2 
                        else
                           FAX(kx,ky,kz) = px * Ekh / mp2
                           FAY(kx,ky,kz) = py * Ekh / mp2
                           FAZ(kx,ky,kz) = pz * Ekh / mp2
                        endif

                        if (isioproc.eq.1) then
                           write (*,*) "Mode"
                           write (*,*) "kappa = ",kx,ky,kz,kappa, sqrt(FAX(kx,ky,kz)**2+FAY(kx,ky,kz)**2+FAZ(kx,ky,kz)**2)
                           write (*,*) "Amplitudes - A"
                           write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                           write (*,*) "Frequencies"
                           write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                           write (*,*) "TAT"
                           write (*,*) TAT(kx,ky,kz), TAP(kx,ky,kz)
                           write (*,*) "Amplitudes - AA"
                           write (*,*) FPXX(kx,ky,kz), FPYX(kx,ky,kz), FPZX(kx,ky,kz)
                           write (*,*) FPXY(kx,ky,kz), FPYY(kx,ky,kz), FPZY(kx,ky,kz)
                           write (*,*) FPXZ(kx,ky,kz), FPYZ(kx,ky,kz), FPZZ(kx,ky,kz)
                        endif
                     endif
                  endif   
               enddo
            enddo
         enddo

c     Now let's break symmetry, have to assume high aspect ratio in z for now
         reduced_mode_count = 0
         do kz = 1, zstep - 1
            kzd = dfloat(kz)
            do ky = mode_start, nymodes
               kyd = dfloat(ky)
               do kx = mode_start, nxmodes
                  kxd = dfloat(kx)
                  
                  kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )

                  if (kappa.le.kappaMax) then
                     call blutilrand(rn)
                     FTX(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                     call blutilrand(rn)
                     FTY(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                     call blutilrand(rn)
                     FTZ(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
c     Translation angles, theta=0..2Pi and phi=0..Pi
                     call blutilrand(rn)
                     TAT(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     TAP(kx,ky,kz) = rn*Pi
c     Phases
                     call blutilrand(rn)
                     FPX(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPY(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPZ(kx,ky,kz) = rn*twicePi
                     if (div_free_force.eq.1) then
                        call blutilrand(rn)
                        FPXX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPXY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPXZ(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYZ(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZZ(kx,ky,kz) = rn*twicePi
                     endif
c     Amplitudes (alpha)
                     call blutilrand(rn)
                     thetaTmp      = rn*twicePi
                     cosThetaTmp   = cos(thetaTmp)
                     sinThetaTmp   = sin(thetaTmp)
                     call blutilrand(rn)
                     phiTmp        = rn*Pi
                     cosPhiTmp     = cos(phiTmp)
                     sinPhiTmp     = sin(phiTmp)

                     px = cosThetaTmp * sinPhiTmp
                     py = sinThetaTmp * sinPhiTmp
                     pz =               cosPhiTmp

                     mp2           = px*px + py*py + pz*pz
                     if (kappa .lt. 0.000001) then
                        if (isioproc .eq. 1) then
                           write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                        endif
                        FAX(kx,ky,kz) = zero
                        FAY(kx,ky,kz) = zero
                        FAZ(kx,ky,kz) = zero
                     else
c     Count modes that contribute
                        reduced_mode_count = reduced_mode_count + 1
c     Set amplitudes
                        if (spectrum_type.eq.1) then
                           Ekh        = one / kappa
                        else if (spectrum_type.eq.2) then
                           Ekh        = one / (kappa*kappa)
                        else
                           Ekh        = one
                        endif
                        if (div_free_force.eq.1) then
                           Ekh = Ekh / kappa
                        endif
                        if (moderate_zero_modes.eq.1) then
                           if (kx.eq.0) Ekh = Ekh / two
                           if (ky.eq.0) Ekh = Ekh / two
                           if (kz.eq.0) Ekh = Ekh / two
                        endif
                        if (force_scale.gt.zero) then
                           FAX(kx,ky,kz) = forcing_epsilon * force_scale * px * Ekh / mp2
                           FAY(kx,ky,kz) = forcing_epsilon * force_scale * py * Ekh / mp2
                           FAZ(kx,ky,kz) = forcing_epsilon * force_scale * pz * Ekh / mp2
                        else
                           FAX(kx,ky,kz) = forcing_epsilon * px * Ekh / mp2
                           FAY(kx,ky,kz) = forcing_epsilon * py * Ekh / mp2
                           FAZ(kx,ky,kz) = forcing_epsilon * pz * Ekh / mp2
                        endif

                        if (isioproc.eq.1) then
                           write (*,*) "Mode"
                           write (*,*) "kappa = ",kx,ky,kz,kappa, sqrt(FAX(kx,ky,kz)**2+FAY(kx,ky,kz)**2+FAZ(kx,ky,kz)**2)
                           write (*,*) "Amplitudes - B"
                           write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                           write (*,*) "Frequencies"
                           write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                           write (*,*) "TAT"
                           write (*,*) TAT(kx,ky,kz), TAP(kx,ky,kz)
                           write (*,*) "Amplitudes - BB"
                           write (*,*) FPXX(kx,ky,kz), FPYX(kx,ky,kz), FPZX(kx,ky,kz)
                           write (*,*) FPXY(kx,ky,kz), FPYY(kx,ky,kz), FPZY(kx,ky,kz)
                           write (*,*) FPXZ(kx,ky,kz), FPYZ(kx,ky,kz), FPZZ(kx,ky,kz)
                        endif
                     endif
                  endif   
               enddo
            enddo
         enddo

         if (isioproc .eq. 1) then
            write(*,*) "mode_count = ",mode_count
            write(*,*) "reduced_mode_count = ",reduced_mode_count
            if (spectrum_type.eq.1) then
               write (*,*) "Spectrum type 1"
            else if (spectrum_type.eq.2) then
               write (*,*) "Spectrum type 2"
            else
               write (*,*) "Spectrum type OTHER"
            endif
         endif
      endif
#endif

      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
#if defined(BL_DO_FLCT)
         write(6,flctin)
#endif
         write(6,control)
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
      
      REAL_T Patm, x, pmf_vals(maxspec+3), sum, Y_pmf(maxspec)
      REAL_T Yfu(maxspec), Yox(maxspec)
      REAL_T Y_FU_H2, Y_FU_O2, Y_OX_H2, Y_OX_O2
      REAL_T X_FU_H2, X_FU_O2, X_FU_N2, Xsum
      REAL_T Xfu(maxspec)
      integer i, zone, n, lo(SDIM), hi(SDIM), iN2, iH2, iO2, getZone

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

      REAL_T X_H2(Nzones),X_O2(Nzones),X_N2(Nzones),T(Nzones)
      REAL_T Xt(maxspec)
c     ZONES:            1                   2         3             4           5     6         7
#ifdef SWIRL
c     phi = 0.8
      data X_H2 /  0.0774907749077d0,   zero, 0.0774907749077d0,   zero,    zero,   zero ,   zero /
      data X_O2  /   0.193726937269d0,  0.21d0, 0.193726937269d0,  0.21d0,  0.21d0, 0.21d0,  0.21d0 /
      data T     /        3.d2,          3.d2,         3.d2,        3.d2,    3.d2,   3.d2 ,   3.d2 /
#else
c     phi = 0.7
c      data X_H2 /  0.0684676292501d0,  zero, 0.0684676292501d0,  zero, zero /
c      data X_O2  /  0.195621797857d0, 0.21d0, 0.195621797857d0, 0.21d0, zero /
c     phi = 1.0
      data X_H2 /  0.095023d0,   zero,  0.095023d0,   zero,   zero /
      data X_O2  /  0.190045d0, 0.21d0,  0.190045d0,  0.21d0,  zero /
      data T     / 3.d2,         3.d2,         3.d2,      3.d2 , zero/
#endif

      integer isioproc
      call bl_pd_is_ioproc(isioproc)

      Patm = pamb / 101325.0d0
      iN2 = -1
      do n = 1,Nspec
         call get_spec_name(name, n)
         if (name .eq. 'N2' ) iN2 = n
      end do

c     For probtype == 3, set boundary conditions as constants
         iN2 = -1
         iH2 = -1
         iO2 = -1
         do n = 1,Nspec
            call get_spec_name(name,n)
            if (name .eq. 'N2' ) iN2  = n
            if (name .eq. 'H2') iH2 = n
            if (name .eq. 'O2' ) iO2  = n
         end do
 
         if(phi_in.ge.zero)then
 
         do n=1,Nspec
            Xfu(n) = zero
         enddo
 
         X_FU_H2 = 2.d0*phi_in
         X_FU_O2 = 1.d0
         X_FU_N2 = .79d0/.21d0
 
         Xsum = X_FU_H2 + X_FU_O2 + X_FU_N2
 
         Xfu(iH2) = X_FU_H2/ Xsum
         Xfu(iO2) = X_FU_O2/ Xsum
         Xfu(iN2) = X_FU_N2/ Xsum

         if (isioproc.eq.1) then
            write(6,*)" in bc routine", Xfu(iH2),Xfu(iO2) ,Xfu(iN2)
         endif
 
 
         call FORT_MOLETOMASS(lo,hi,
     &        Xfu, ARLIM(lo),ARLIM(hi),
     &        Yfu,         ARLIM(lo),ARLIM(hi))
 
         else
c     phi = 1.0
         Y_FU_H2 = 0.0551670662321d0
         Y_FU_O2  = 0.220068142419d0
c     phi = 0.7
c         Y_FU_H2 = 0.0392668168889d0
c         Y_FU_O2  = 0.223771589041d0
 
 
         do n=1,Nspec
            Yfu(n) = zero
         enddo
         Yfu(iH2) = Y_FU_H2
         Yfu(iO2)  = Y_FU_O2
         Yfu(iN2)  = one - Y_FU_H2 - Y_FU_O2
 
         endif
 
         zone = 1
         do n=1,Nspec
            Y_bc(n-1,zone) = Yfu(n)
         end do
         T_bc(zone) = T_in
         u_bc(zone) = zero
         v_bc(zone) = V_in
 
         call FORT_RHOfromPTY(lo,hi,
     &                        rho_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        T_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        Y_bc(0,zone),ARLIM(lo),ARLIM(hi), Patm)
         call FORT_HMIXfromTY(lo,hi,
     &                        h_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        T_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        Y_bc(0,zone),ARLIM(lo),ARLIM(hi))
 

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
c     getZone = BL_VOLUME
c     if ( ((direction.eq."z") .and. (z.le.domnlo(3))).or.
c    &     ((direction.eq."y") .and. (y.le.domnlo(2))).or.
c    &     ((direction.eq."x") .and. (x.le.domnlo(1)))) then


            getZone = 1

c      endif

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
      REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h, dx(SDIM)
      logical getuvw

#include "cdwrk.H"
#include "htdata.H"
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

         if (getuvw) then
         if (probtype.eq.3) then
            u = zero
            v = zero
            w =  V_in + (time-tbase_control)*dV_control
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
#include"probdata.H"
      turbSclT = 0.d0
      eta = r/Ro
      if(r.lt.Ro)then
#ifdef SWIRL
         turbSclT = V_in*(0.10278d0 + 0.34958d0 *eta - 0.86045d0 *eta**2
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
         turbSclR=V_in*(0.086871d0 + 1.5445d0*eta - 17.316d0*eta**2 + 79.469d0*eta**3
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
c         turbSclZ = V_in*(-0.020318d0 + 2.18d0*eta - 18.157d0*eta**2 + 66.538d0*eta**3 - 
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
      convVelSAVE = convVel
      convVel = 15.d0
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
#endif
      enddo            

         do k= lo(3),hi(3)
            z = (float(k)+.5)*delta(3)+domnlo(3)            
            do j = lo(2), hi(2)
               y = (float(j)+5)*delta(2)+domnlo(2)
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
      convVel = convVelSAVE
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
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),
     &                 ftrac(i,j,k,1).gt.flametracval)
               enddo
            enddo
         enddo
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
c                 ax = abs(temperature(i+1,j,k,1) - temperature(i,j,k,1))
c                 ay = abs(temperature(i,j+1,k,1) - temperature(i,j,k,1))
c                 az = abs(temperature(i,j,k+1,1) - temperature(i,j,k,1))
c                 ax = MAX(ax,abs(temperature(i,j,k,1) - temperature(i-1,j,k,1)))
c                 ay = MAX(ay,abs(temperature(i,j,k,1) - temperature(i,j-1,k,1)))
c                 az = MAX(az,abs(temperature(i,j,k,1) - temperature(i,j,k-1,1)))
c                 aerr = max(ax,ay,az)
c                 tag(i,j,k) = cvmgt(set,tag(i,j,k),aerr.ge.tempgrad)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),temperature(i,j,k,1).lt.temperr)
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
     &                         ABS(vort(i,j,k,1)).ge.vorterr*2.d0**level)
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
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
      REAL_T  cs,sg,r,scaler,scalet
      REAL_T  turbSclX,turbSclY,turbSclZ,turbSclR,turbSclT

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
         dt_flct = time - tbase_control
         t_flct = zbase_control + V_in*dt_flct + dV_control*dt_flct**2
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
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)

                  call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(1,1,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 1) then
c                       vel(i,j,k,1) = u + uflct(1,j,k)*turbSclX(x,y,z)*turb_scale
                        vel(i,j,k,1) = u + uflct(1,j,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(1,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 1) then
c                       vel(i,j,k,2) = v + vflct(1,j,k)*turbSclY(x,y,z)*turb_scale
                        vel(i,j,k,2) = v + vflct(1,j,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(1,1,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 1) then
c                       vel(i,j,k,3) = w + wflct(1,j,k)*turbSclZ(x,y,z)*turb_scale
                        vel(i,j,k,3) = w + wflct(1,j,k)*turb_scale
                     endif
#endif
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
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do i = lo(1), hi(1)
                  x = (float(i)+.5)*delta(1)+domnlo(1)

                  call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(2,1,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 2) then
c                       vel(i,j,k,1) = u + uflct(i,1,k)*turbSclX(x,y,z)*turb_scale
                        vel(i,j,k,1) = u + uflct(i,1,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 2) then
c                       vel(i,j,k,2) = v + vflct(i,1,k)*turbSclY(x,y,z)*turb_scale
                        vel(i,j,k,2) = v + vflct(i,1,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,1,3).eq.EXT_DIR) then
                        vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 2) then
c                       vel(i,j,k,3) = w + wflct(i,1,k)*turbSclZ(x,y,z)*turb_scale
                        vel(i,j,k,3) = w + wflct(i,1,k)*turb_scale
                     endif

#endif
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

                  call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(2,2,1).eq.EXT_DIR) then
                        vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 2) then
c                       vel(i,j,k,1) = u + uflct(i,1,k)*turbSclX(x,y,z)*turb_scale
                        vel(i,j,k,1) = u + uflct(i,1,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 2) then
c                       vel(i,j,k,2) = v + vflct(i,1,k)*turbSclY(x,y,z)*turb_scale
                        vel(i,j,k,2) = v + vflct(i,1,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 2) then
c                       vel(i,j,k,3) = w + wflct(i,1,k)*turbSclZ(x,y,z)*turb_scale
                        vel(i,j,k,3) = w + wflct(i,1,k)*turb_scale
                     endif
#endif
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

                  call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(3,1,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,1) = u + uflct(i,j,1)*turb_scale
                     endif
#endif
                  endif

                  if (bc(3,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 3) then
                          vel(i,j,k,2) = v + vflct(i,j,1)*turb_scale
                     endif
#endif
                  endif

                  if (bc(3,1,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,3) = w + wflct(i,j,1)*turb_scale
                     endif
#endif
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

                  call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                  if (bc(3,2,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 3) then
c                       vel(i,j,k,1) = u + uflct(i,j,1)*turbSclX(x,y,z)*turb_scale
                        vel(i,j,k,1) = u + uflct(i,j,1)*turb_scale
                     endif
#endif
                  endif

                  if (bc(3,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 3) then
c                       vel(i,j,k,2) = v + vflct(i,j,1)*turbSclY(x,y,z)*turb_scale
                        vel(i,j,k,2) = v + vflct(i,j,1)*turb_scale
                     endif
#endif
                  endif

                  if (bc(3,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 3) then
c                       vel(i,j,k,3) = w + wflct(i,j,1)*turbSclZ(x,y,z)*turb_scale
                        vel(i,j,k,3) = w + wflct(i,j,1)*turb_scale
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
#endif      

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
                     xvel(i,j,k) = u + uflct(1,j,k)*turb_scale
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
                     xvel(i,j,k) = u + uflct(1,j,k)*turb_scale
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
                     xvel(i,j,k) = u + uflct(i,1,k)*turb_scale
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
                     xvel(i,j,k) = u + uflct(i,1,k)*turb_scale
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
                     xvel(i,j,k) = u + turb_scale
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
                     xvel(i,j,k) = u + uflct(i,j,1)*turb_scale
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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#if defined(BL_DO_FLCT)
#include "INFL_FORCE_F.H"
#endif
      
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
                     yvel(i,j,k) = v + vflct(1,j,k)*turb_scale
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
                     yvel(i,j,k) = v + vflct(1,j,k)*turb_scale
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
                     yvel(i,j,k) = v + vflct(i,1,k)*turb_scale
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
                     yvel(i,j,k) = v + vflct(i,1,k)*turb_scale
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
                     yvel(i,j,k) = v + vflct(i,j,1)*turb_scale
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
                     yvel(i,j,k) = v + vflct(i,j,1)*turb_scale
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

#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#if defined(BL_DO_FLCT)
#include "INFL_FORCE_F.H"
#endif
      
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
                     zvel(i,j,k) = w + wflct(1,j,k)*turb_scale
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
                     zvel(i,j,k) = w + wflct(1,j,k)*turb_scale
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
                     zvel(i,j,k) = w + wflct(i,1,k)*turb_scale
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
                     zvel(i,j,k) = w + wflct(i,1,k)*turb_scale
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
                     zvel(i,j,k) = w + wflct(i,j,1)*turb_scale
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
                     zvel(i,j,k) = w + wflct(i,j,1)*turb_scale
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

c     if (add_turb .and. zturb_top .gt. zero) then

c        add_turb_var = 1
c        turb_scale_var = turb_scale

c        turb_length(1) = domnhi(1)-domnlo(1)
c        turb_length(2) = domnhi(2)-domnlo(2)
c        turb_length(3) = zturb_top - zturb_bot
c        
c        box_offset(1) = 0
c        box_offset(2) = 0
c        box_offset(3) = int(zturb_bot/delta(3))
c     endif

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
      
c     Blend goes from 1 to 0 in r at wire_r, over a width of wire_w  
c     blend(r) = one - half*(one+TANH((r-wire_r)/(fourth*wire_w)))
      return

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
               rad(i,j,k) = half*(one-tanh(two*(r-1.5*hotspot_w)/hotspot_w))
     &                                       *hotspot_htc*(hotspot_T - T(i,j,k))

#else
c              r = SQRT((xx-wire_x)**2 + (zz-wire_z)**2 + MAX(zero,ABS(yy)-wire_y)**2)
c              rad(i,j,k) = -blend(r)*wire_h*MAX(zero,(wire_T - T(i,j,k)))
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

      REAL_T     Lx, Ly, Lz, Lmin, HLx, HLy, HLz
      REAL_T     kappa, kappaMax

#include "probdata.H"
#include "cdwrk.H"
#include "bc.H"

#ifdef DO_LMC_FORCE
#include "forcedata.H"
#endif

c
c     ::::: local variables
c
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer a2, a3, a4, a5
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f1, f2, f3
      REAL_T  twicePi
      REAL_T  force_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc
      integer nXvel, nYvel, nZvel, nRho, nTrac

      if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
         write(*,*) "pseudo_gravity::dV_control = ",dV_control
      endif

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      ilo = istate_l1
      jlo = istate_l2
      klo = istate_l3
      ihi = istate_h1
      jhi = istate_h2
      khi = istate_h3

c     Assumes components are in the following order
      nXvel = 1
      nYvel = 2
      nZvel = 3
      nRho  = 4
      nTrac = 5

      if (scomp.eq.0) then
c     Do velocity forcing
         if (probtype.eq.3) then
#ifdef DO_LMC_FORCE
c     Homogeneous Isotropic Turbulence
            twicePi=two*Pi
            
c     Adjust z offset for probtype 15
            if (time_offset.gt.(-half)) then
               force_time = time + time_offset
            else
               force_time = time
            endif

            Lx = domnhi(1)-domnlo(1)
            Ly = domnhi(2)-domnlo(2)
            Lz = domnhi(3)-domnlo(3)

            if (hack_lz.eq.1) then 
               Lz = Lz/two
            endif
         
            Lmin = min(Lx,Ly,Lz)
            kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
            nxmodes = nmodes*int(0.5+Lx/Lmin)
            nymodes = nmodes*int(0.5+Ly/Lmin)
            nzmodes = nmodes*int(0.5+Lz/Lmin)

            xstep = int(Lx/Lmin+0.5)
            ystep = int(Ly/Lmin+0.5)
            zstep = int(Lz/Lmin+0.5)

            if (forcing_twice_wavelength.eq.1) then
               HLx = Lx/two
               HLy = Ly/two
               HLz = Lz/two
            else
               HLx = Lx
               HLy = Ly
               HLz = Lz
            endif

            do k = klo, khi
               z = xlo(3) + hz*(float(k-klo) + half)
               do j = jlo, jhi
                  y = xlo(2) + hy*(float(j-jlo) + half)
                  do i = ilo, ihi
                     x = xlo(1) + hx*(float(i-ilo) + half)
                     f1 = zero
                     f2 = zero
                     f3 = zero
                     do kz = mode_start*zstep, nzmodes, zstep
                        kzd = dfloat(kz)
                        do ky = mode_start*ystep, nymodes, ystep
                           kyd = dfloat(ky)
                           do kx = mode_start*xstep, nxmodes, xstep
                              kxd = dfloat(kx)
                              kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                              if (kappa.le.kappaMax) then
                                 xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                                 if (div_free_force.eq.1) then
                                    f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) 
     &                                   -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                                    f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) 
     &                                   -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                                    f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) 
     &                                   -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                                 else
                                    f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                    f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                    f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                 endif
                              endif
                           enddo
                        enddo
                     enddo
                     do kz = 1, zstep - 1
                        kzd = dfloat(kz)
                        do ky = mode_start, nymodes
                           kyd = dfloat(ky)
                           do kx = mode_start, nxmodes
                              kxd = dfloat(kx)
                              kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                              if (kappa.le.kappaMax) then
                                 xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                                 if (div_free_force.eq.1) then
                                    f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) 
     &                                   -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                                    f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) 
     &                                   -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                                    f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) 
     &                                   -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                                 else
                                    f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                    f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                    f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                 endif
                              endif
                           enddo
                        enddo
                     enddo
                     if (use_rho_in_forcing.eq.1) then
                        force(i,j,k,nXvel) = f1*rho(i,j,k)
                        force(i,j,k,nYvel) = f2*rho(i,j,k)
                        force(i,j,k,nZvel) = f3*rho(i,j,k)
                     else
                        force(i,j,k,nXvel) = f1
                        force(i,j,k,nYvel) = f2
                        force(i,j,k,nZvel) = f3
                     endif
                  enddo
               enddo
            enddo
#endif
c     Default to gravity...
         elseif (abs(gravity).gt.0.0001) then
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = gravity*rho(i,j,k)
                  enddo
               enddo
            enddo
c     else to zero
         else
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = zero
                  enddo
               enddo
            enddo
         endif
c     Let's add the pseudo gravity afterwards...
         if (pseudo_gravity.eq.1) then
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nZvel) = force(i,j,k,nZvel)  + dV_control*rho(i,j,k)
                  enddo
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
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else if (n.eq.nTrac) then
c     Tracer
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else
c     Other scalar
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            endif
         enddo
      endif

      end

