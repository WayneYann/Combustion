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

#define SDIM 2

#if defined(BL_T3E) || defined(BL_AIX) || defined(BL_IRIX64) || defined(BL_CRAYX1)
#define OPEN_APPEND position='append'
#else
#define OPEN_APPEND access='append'
#endif

      subroutine FORT_ACTIVECONTROL(coft,time,dt,myproc,step,restart)
      implicit none
      REAL_T coft,time,dt
      integer myproc,step,restart
      call bl_abort('Should not in FORT_ACTIVECONTROL...')
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
#include "visc.H"
#include "conp.H"

      integer i, n
      character*(maxspnml) spName
      integer lo(SDIM),hi(SDIM)
      data lo /1,1/
      data hi /1,1/

      namelist /fortin/ vorterr, temperr, adverr, tempgrad, flametracval,
     &                  probtype, max_temp_lev, max_vort_lev, 
     &                  max_trac_lev,T_amb,
     &                  v_strength,v_width,v_xcen,v_ycen,v_cl_x,
     &                  phi_in,r_hot,x_hot,y_hot,hot_wid,T_hot,X_nitrogen
      namelist /heattransin/ pamb, dpdt_factor

c      Build `probin' filename -- the name of file containing fortin namelist.
c
      integer maxlen
      parameter (maxlen=256)

      character probin*(maxlen)

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

c     Set default values for namelisted variables
      vorterr = 1.d20
      temperr = 1.d0
      adverr =  1.d0
      tempgrad = 1.d20
      flametracval = 1.d20
      probtype = 0
      max_temp_lev = 100
      max_vort_lev = 100
      max_trac_lev = 100
      v_strength = 0.d0
      v_width = 0.d0
      v_xcen = 0.d0
      v_ycen = 0.d0
      v_cl_x = 0.d0
      phi_in = 0.d0
      r_hot = 0.d0
      x_hot = 0.d0
      y_hot = 0.d0
      hot_wid = 1.d0
      T_hot = 1.d3
      X_nitrogen = 0.d0
      pamb = 101325.d0
      dpdt_factor = 0.3d0
      
      read(untin,fortin)
c      write(6,fortin)

      read(untin,heattransin)
c      write(6,heattransin)

   10 close(unit=untin)

c     Load domain dimensions into common, and set up boundary functions
      domnlo(1) = problo(1)
      domnlo(2) = problo(2)
      domnhi(1) = probhi(1)
      domnhi(2) = probhi(2)

      iN2 = -1
      iH2 = 1
      iO2 = -1
      do n = 1,Nspec
         call get_spec_name(spName,n)
         if (spName .eq. 'N2' ) iN2 = n
         if (spName .eq. 'H2' ) iH2 = n
         if (spName .eq. 'O2' ) iO2 = n
      end do

      do n=1,Nspec
         X_comp(n) = zero
      enddo
      X_comp(iN2) = X_nitrogen
      X_comp(iH2) = phi_in/(phi_in+0.5d0)*(1.d0-X_nitrogen)
      X_comp(iO2) = 1.d0 - X_nitrogen -  X_comp(iH2)
      call FORT_MOLETOMASS(lo,hi,
     &     X_comp, ARLIM(lo),ARLIM(hi),
     &     Y_comp, ARLIM(lo),ARLIM(hi))

      call setupbc()
      bcinit = .true.

      end

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

      subroutine setupbc()
      implicit none
#include "probdata.H"
#include "htdata.H"
      integer zone, n
      integer lo(SDIM),hi(SDIM)
      data lo /1,1/
      data hi /1,1/
      REAL_T Patm

      zone = 0
      do n=1,Nspec
         Y_bc(n-1,zone) = Y_comp(n)
      end do
      T_bc(zone) = T_amb
      u_bc(zone) = zero
      v_bc(zone) = zero

      Patm = pamb / 101325.d0
      
      call FORT_RHOfromPTY(lo,hi,
     &                     rho_bc(zone),ARLIM(lo),ARLIM(hi),
     &                     T_bc(zone),ARLIM(lo),ARLIM(hi),
     &                     Y_bc(0,zone),ARLIM(lo),ARLIM(hi), Patm)
      call FORT_HMIXfromTY(lo,hi,
     &                     h_bc(zone),ARLIM(lo),ARLIM(hi),
     &                     T_bc(zone),ARLIM(lo),ARLIM(hi),
     &                     Y_bc(0,zone),ARLIM(lo),ARLIM(hi))

      end

         
c ::: -----------------------------------------------------------
      
      integer function getZone(x, y)
      implicit none
#include "probdata.H"
      REAL_T x, y
      getZone = 1
      end

c ::: -----------------------------------------------------------
      
      subroutine bcfunction(orient, x, y, time, u, v, rho, Yl, T, h, dx)
      implicit none
      integer orient
      REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(SDIM)
#include "htdata.H"
#include "probdata.H"

      integer n, getZone, zone
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
      else
          print*, 'No boundary condition for orientation = ', orient
          call bl_abort(' ')
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
#include "probdata.H"

      integer i, j, n
      REAL_T  x, y
      integer slo(SDIM), shi(SDIM)
      character*(maxspnml) name
      REAL_T Patm
      integer spec, getZone, zone
      
      integer zL, zR, zH
      REAL_T eta, xi, xblend, mdotb, mdotL, mdotR, r1sq, r2sq, ylo, yhi
      REAL_T pmf_H_vals(maxspec+3), sum, Y_pmf(maxspec), rblob, theta

      if ((Temp.gt.0).neqv.(RhoH.gt.0)) then
         call bl_abort('Need both Temp and RhoH, or neither')
      end if
      
      if ((Temp .LT. 0) .OR. (RhoH .LT. 0)) then
         call bl_abort('No IC''s for system without T, RhoH')
      endif

c
c     Set vorticity into pressure field
c
      do j = ARG_L2(press),ARG_H2(press)
         y = float(j)*delta(2)+domnlo(2)
         do i = ARG_L1(press),ARG_H1(press)
            x = float(i)*delta(1)+domnlo(1)
            r1sq = (x-v_xcen)**2              + (y-v_ycen)**2
            r2sq = (x-(two*v_cl_x-v_xcen))**2 + (y-v_ycen)**2
            press(i,j) = 2.d0*v_strength/((1.d0-dexp(-1.d0))*v_width)*
     &           (dexp(-r1sq/v_width**2)-dexp(-r2sq/v_width**2)) 
         end do
      end do
      
      do j = lo(2), hi(2)
         y = (float(j)+.5)*delta(2)+domnlo(2)         
         do i = lo(1), hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            rblob = sqrt( (x_hot-x)**2 + (y_hot-y)**2 )

            scal(i,j,Temp) = 300.d0+ (T_hot-300.d0)*dexp(-(rblob/hot_wid)**2)
            do n = 1,Nspec
               scal(i,j,firstSpec+n-1) = Y_comp(n)
            end do
            
            vel(i,j,1) = zero
            vel(i,j,2) = zero

            scal(i,j,Trac) = fourth*(press(i,j  ) + press(i+1,j)
     &           +                   press(i,j+1) + press(i+1,j+1))

         enddo
      enddo
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
      call bl_abort('DENERROR: should no be here')
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
               tag(i,j) = cvmgt(set,tag(i,j),ftrac(i,j,1).gt.flametracval)
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
               ax = abs(temperature(i+1,j,1) - temperature(i-1,j,1))
               ay = abs(temperature(i,j+1,1) - temperature(i,j-1,1))
               aerr = max(ax,ay)
               tag(i,j) = cvmgt(set,tag(i,j),aerr.ge.tempgrad)
c              tag(i,j) = cvmgt(set,tag(i,j),temperature(i,j,1).le.tempgrad)
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
               tag(i,j) = cvmgt(set,tag(i,j),ABS(v).ge.vorterr)
            end do
         end do
      endif

      if (level .lt. max_vort_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j) = cvmgt(set,tag(i,j),ABS(vort(i,j,1)).ge.vorterr*2.d0**level)
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
           do i = lo(1), domlo(1)-1
             do j = lo(2), hi(2)
               adv(i,j) = zero
             enddo
           enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
           do i = domhi(1)+1, hi(1)
             do j = lo(2), hi(2)
               adv(i,j) = zero
             enddo
           enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
           do j = lo(2), domlo(2)-1
              do i = lo(1), hi(1)
                 adv(i,j) = zero
              enddo
           enddo
      endif            

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
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

#include "probdata.H"

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

#include "probdata.H"
      
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
     &                          xlo,time,bc,id )
      implicit none
      integer DIMDEC(rhoY), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM), id
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoY(DIMV(rhoY))

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
           divu(i,j) = cvmgt(one,zero,x.GE.0.5)
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
      end
      
      subroutine FORT_RADLOSS(lo,hi,rad,DIMS(rad),
     &                        T,DIMS(T),Y,DIMS(Y),dx,Patm,time)
      implicit none
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
