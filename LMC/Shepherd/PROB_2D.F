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

      namelist /fortin/ flowrate_fu_CH4, flowrate_fu_N2,
     &     flowrate_fu_O2,T_fu, flowrate_fu_H2, flowrate_fu_NO, R_fu,
     &     flowrate_ox_O2, flowrate_ox_N2, T_ox, R_ox_1, R_ox_2,
     &     x_hot, y_hot_lo, y_hot_hi, T_lip, bgBL, flBL, ystandoff,
     &     T_hot, temperr, adverr, vorterr, tempgrad, max_temp_lev,
     &     max_vort_lev, flametracval, y_wall2, y_wall3, T_wall1,
     &     T_wall2, T_wall3, Trad_bg, time_rad, dt_rad,
     &     refine_nozzle, nozzlerefwidth, nozzlerefheight, 
     &     nozzlerefloc, probtype, max_nozzle_lev
      namelist /heattransin/ pamb, dpdt_factor
c
c      Build 'probin' filename -- the name of file containing fortin namelist.
c
      integer maxlen, i
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

c     Set defaults
      flowrate_fu_CH4 = one
      flowrate_fu_N2 = one
      flowrate_fu_O2 = zero
      flowrate_fu_H2 = zero
      flowrate_fu_NO = zero
      T_fu = 300d0
      R_fu = .01d0
      flowrate_ox_O2 = one
      flowrate_ox_N2 = one
      T_ox = 300d0
      R_ox_1 = .02d0
      R_ox_2 = .03d0
      x_hot = R_fu
      y_hot_hi = R_fu
      y_hot_lo = zero
      T_hot = 2.d3
      T_lip = 1.5d3
      T_wall1 = 300.d0
      T_wall2 = 300.d0
      T_wall3 = 300.d0
      y_wall2 = zero
      y_wall3 = zero
      Trad_bg = 300.d0
      time_rad = zero
      dt_rad = 1.e-9
      bgBL = one
      flBL = one

      temperr = zero
      adverr = 1.e20
      vorterr = 5.e3
      tempgrad  = 50.0d0
      max_temp_lev = 0
      max_vort_lev = 0
      flametracval = .005
      
      refine_nozzle = 0
      nozzlerefwidth = zero
      nozzlerefheight = zero
      nozzlerefloc = zero
      max_nozzle_lev = 0
      probtype = 0

      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      
      read(untin,fortin)
      
      pamb = 101325.0
      dpdt_factor = 0.3
      read(untin,heattransin)

   10 close(unit=untin)

c     Load domain dimensions into common, and set up boundary functions
      domnlo(1) = problo(1)
      domnlo(2) = problo(2)
      domnhi(1) = probhi(1)
      domnhi(2) = probhi(2)

c     Precompute some messy coeffs
      a_ox = (R_ox_2**2 - R_ox_1**2)/(4*log(R_ox_1/R_ox_2))
      b_ox = (log(R_ox_2)*R_ox_1**2 - log(R_ox_1)*R_ox_2**2) /
     &     (4*log(R_ox_1/R_ox_2))
      c_ox = (8.d0/(Pi*6.d7)) /
     &     (R_ox_2**2*(-4*a_ox+8*b_ox+R_ox_2**2+8*a_ox*log(R_ox_2))
     &     -R_ox_1**2*(-4*a_ox+8*b_ox+R_ox_1**2+8*a_ox*log(R_ox_1)))
      
      call setupbc()
      bcinit = .true.
      
      end

c ::: -----------------------------------------------------------

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
      
      subroutine setupbc()
      implicit none
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
#include "htdata.H"
      
      REAL_T Patm, X(0:maxspec-1), Y(0:maxspec-1), tot
      REAL_T flowrate_fu_tot, flowrate_ox_tot
      integer n, zone, iN2, lo(SDIM), hi(SDIM)
      character*(maxspnml) name
      data lo / 1, 1 /
      data hi / 1, 1 /
      
      Patm = pamb / 101325.0
      iN2 = -1
      do n = 0,Nspec-1
         call get_spec_name(name, n+1)
         if (name .eq. 'N2' ) iN2 = n
      end do
c
c	Glarborg and Bessler cases:
c
c	Boundary is divided into 5 zones
      if (probtype.eq.0) then	      
c     Compute inlet mass fracs based on input flowrates
c
c     ZONE 1: Fuel stream
c     Needs: flowrate_fu_CH4, flowrate_fu_N2, flowrate_fu_O2, flowrate_fu_H2,
c     T_fu, R_fu (flowrates in ml/min, R_fu in mm)
c     
c     Note: 6.e7 = 60 s/min * 1.e6 ml/m^3
         zone = 1
         tot = zero
         flowrate_fu_tot = flowrate_fu_CH4 + flowrate_fu_N2
     &        + flowrate_fu_O2 + flowrate_fu_H2 + flowrate_fu_NO
         do n = 0,Nspec-1
            call get_spec_name(name, n+1)
            if (name .EQ. 'CH4') then
               X(n) = flowrate_fu_CH4/flowrate_fu_tot
            else if (name .EQ. 'O2') then
               X(n) = flowrate_fu_O2/flowrate_fu_tot
            else if (name .EQ. 'H2') then
               X(n) = flowrate_fu_H2/flowrate_fu_tot
            else if (name .EQ. 'NO') then
               X(n) = flowrate_fu_NO/flowrate_fu_tot
            else
               X(n) = zero
            endif
            tot = tot + X(n)
         end do
         X(iN2) = one - tot

         CALL CKXTY (X, IWRK, RWRK, Y_bc(0,zone))

         T_bc(zone) = T_fu
         u_bc(zone) = zero
         v_bc(zone) = (flowrate_fu_tot/6.e7) / (Pi*R_fu**2)
c         
c     ZONE 2: Oxidizer stream
c     Needs: flowrate_ox_O2, flowrate_ox_N2, T_ox, R_ox_1, R_ox_2
c     (flowrate in ml/min, Rs in mm)
c
         zone = 2
         tot = zero
         flowrate_ox_tot = flowrate_ox_O2 + flowrate_ox_N2
         do n = 0,Nspec-1
            call get_spec_name(name, n+1)
            if (name .EQ. 'O2') then
               X(n) = flowrate_ox_O2/flowrate_ox_tot
            else
               X(n) = zero
            endif
            tot = tot + X(n)
         end do
         X(iN2) = one - tot

         CALL CKXTY (X, IWRK, RWRK, Y_bc(0,zone))

         T_bc(zone) = T_ox
         u_bc(zone) = zero
         v_bc(zone) = (flowrate_ox_tot/6.e7) / (Pi*(R_ox_2**2 - R_ox_1**2))
c         
c     ZONE 3: Dead zone, like 2 but no inflow, different T
c
         zone = 3
         do n = 0,Nspec-1
            Y_bc(n,zone) = Y_bc(n,2)
         end do
         T_bc(zone) = T_lip
         u_bc(zone) = u_bc(2)
         v_bc(zone) = zero
c         
c     ZONE 4: Hot spot, like 2 but hotter
c
         zone = 4
         do n = 0,Nspec-1
            Y_bc(n,zone) = Y_bc(n,2)
         end do
         T_bc(zone) = T_hot
         u_bc(zone) = u_bc(2)
         v_bc(zone) = v_bc(2)
c
c     ZONE 5: Outter wall
c
         zone = 5
         do n = 0,Nspec-1
            Y_bc(n,zone) = Y_bc(n,2)
         end do
         T_bc(zone) = T_wall3
         u_bc(zone) = zero
         v_bc(zone) = zero
      else
         write(*,*) 'PROB_2D::setupbc: Dont know about probtype = ',
     &        probtype
         stop
      end if
c
c     Get rho and h for all zones
c
      do zone=1,Nzones
         call FORT_RHOfromPTY(lo,hi,
     &                        rho_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        T_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        Y_bc(0,zone),ARLIM(lo),ARLIM(hi), Patm)
         call FORT_HMIXfromTY(lo,hi,
     &                        h_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        T_bc(zone),ARLIM(lo),ARLIM(hi),
     &                        Y_bc(0,zone),ARLIM(lo),ARLIM(hi))
      end do
      end
         
c ::: -----------------------------------------------------------
      
      integer function getZone(x, y)
      implicit none
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      REAL_T x, y, epsilon

      epsilon = 1.d-9
      getZone = 2
      if (y .le. domnlo(2)+epsilon) then
         if (x.le.R_fu) then
            getZone = 1
         else if ((x.gt.R_fu).and.(x.le.R_ox_1)) then
            getZone = 3
         else if ((x.gt.R_ox_1).and.(x.le.R_ox_2)) then
            getZone = 2
         endif
      else if (x.ge.domnhi(1)-epsilon) then
         getZone = 5
      else
         if ( (y .le. y_hot_hi) .AND.
     &        (y .ge. y_hot_lo) .AND.
     &        (x .le. x_hot)    ) getZone = 4
      end if
      end
      
c ::: -----------------------------------------------------------
      
      subroutine bcfunction(orient, x, y, time, u, v, rho, Yl, T, h, dx)
      implicit none
      integer orient
      REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(SDIM)
#include "htdata.H"
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"
      REAL_T flowrate_fu_tot, flowrate_ox_tot, Patm, xleft, xrght
	REAL_T uavFU,R1FU,R2FU,FlowFU
	REAL_T uavOX,R1OX,R2OX,FlowOX, OXInt, OXIntR1, OXIntR2
      integer n, getZone, zone, lo(SDIM), hi(SDIM)
      data lo / 1, 1 /
      data hi / 1, 1 /
      
c     Compute the average velocity in the fuel stream between R1FU and R2FU
c     Note: If R2FU > R_fu, the average is computed assuming the remainder
c     extends into a no-flow zone (ie, the fuel tube is wider than dx(1))
      uavFU(R1FU,R2FU,FlowFU) = 
     &     (R_fu**2 * (MIN(R_fu,R2FU)**2 - R1FU**2)
     &     -   half * (MIN(R_fu,R2FU)**4 - R1FU**4) )
     &     *FlowFU/(Pi * BL_REAL_E(3.0,7) * (R1FU+R2FU) * dx(1) * R_fu**4)
c     uavFU(R1FU,R2FU,FlowFU) = FlowFU/(BL_REAL_E(3.0,7)*Pi*R_fu**4)
c     &	*(R_fu**2 - (half*(R1FU+R2FU))**2)

c     Compute the average velocity in the oxidizer stream between R1OX and R2OX
c     Note: If R1OX < R_ox_1, or if R2OX > R_ox_2, the average is computed assuming
c            the remainder extends into a no-flow zone (ie, the fuel tube is
c            wider than dx(1))
      OXInt(OXIntR1,OXIntR2) = BL_REAL_E(16.0,-1)*
     &     (OXIntR2**2*(OXIntR2**2+eight*b_ox+four*a_ox*(two*LOG(OXIntR2)-one))
     &     -OXIntR1**2*(OXIntR1**2+eight*b_ox+four*a_ox*(two*LOG(OXIntR1)-one)))

      uavOX(R1OX,R2OX,FlowOX) =
     &     OXInt(MAX(R_ox_1,R1OX),MIN(R_ox_2,R2OX))/OxInt(R_ox_1,R_ox_2)
     &     *(FlowOX/(BL_REAL_E(6.0,7)*Pi*(R1OX+R2OX)*dx(1))) 

      if (.not. bcinit) then
         write(6,*) 'Need to initialize boundary condition function'
         stop
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
c     Overwrite constant inflow velocities with profiles
         flowrate_fu_tot = flowrate_fu_CH4 + flowrate_fu_N2
     &        + flowrate_fu_O2 + flowrate_fu_H2 + flowrate_fu_NO
         flowrate_ox_tot = flowrate_ox_O2 + flowrate_ox_N2

         if ((probtype .eq. 0)) then
            if (zone.eq.1) then
               v = flowrate_fu_tot*(R_fu**2 - x**2)/(3.d7*Pi*R_fu**4)
            else if (zone.eq.2) then
               v = flowrate_ox_tot
     &              * (fourth*c_ox*x**2 + a_ox*c_ox*log(x) + b_ox*c_ox)
            endif
         else
            write(*,*) 'PROB_2D::bcfunction: Dont know about probtype = ',
     &           probtype
            stop
         end if
c     Overwrite constant T on wall with T profile, then recalc rho and h
         if (zone.eq.5) then
            Patm = pamb / 101325.0
            if (y.lt.y_wall2) then
               T = T_wall1
     &              + (T_wall2-T_wall1)*(y-domnlo(2))/(y_wall2-domnlo(2))
            else if (y.lt.y_wall3) then
               T = T_wall2
     &              + (T_wall3-T_wall2)*(y-y_wall2)/(y_wall3-y_wall2)
            else
               T = T_wall3
            end if
            call FORT_RHOfromPTY(lo,hi,
     &                           rho,ARLIM(lo),ARLIM(hi),
     &                           T,  ARLIM(lo),ARLIM(hi),
     &                           Yl, ARLIM(lo),ARLIM(hi), Patm)
            call FORT_HMIXfromTY(lo,hi,
     &                           h, ARLIM(lo),ARLIM(hi),
     &                           T, ARLIM(lo),ARLIM(hi),
     &                           Yl,ARLIM(lo),ARLIM(hi))
         end if 
      else
         write(6,*) 'No boundary condition for orientation = ', orient
         stop
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
     &	 	               vel,scal,DIMS(scal),press,DIMS(press),
     &                         delta,xlo,xhi)
      implicit none
      integer    level, nscal, nvel
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(scal)
      integer    DIMDEC(press)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     time, delta(SDIM)
      REAL_T     vel(DIMV(scal),SDIM)
      REAL_T    scal(DIMV(scal),nscal)
      REAL_T   press(DIMV(press))

#include "htdata.H"
#include "cdwrk.H"
#include "bc.H"
#include "probdata.H"

c     ::::: local variables
      integer i, j, n, ii, jj
      REAL_T  x, y, y1, y2
      REAL_T  u, v, rho, Xl(maxspec), Yl(maxspec), T, h, Patm
      REAL_T xcrit
      REAL_T eta, eta0, xi, sum
      integer getZone, zone, zH, zL, zR, iO2, iN2
      character*(maxspnml) name
      
      if ((Temp.gt.0).neqv.(RhoH.gt.0)) then
         write(6,*) 'Need both Temp and RhoH, or neither'
         stop
      end if
      
      if ((Temp .LT. 0) .OR. (RhoH .LT. 0)) then
         write(6,*) 'No IC''s for system without T, RhoH'
         stop
      endif

      do n = 1,Nspec
         call get_spec_name(name, n)
         if (name .eq. 'N2' ) iN2 = n
         if (name .eq. 'O2' ) iO2 = n
      end do

      if ((probtype.eq.0)) then
         do j = lo(2), hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)

               scal(i,j,Trac) = zero
               scal(i,j,Temp) = T_fu
               do n = 1,Nspec
                  Xl(n) = zero
               end do 
               Xl(iN2) = 0.79d0
               Xl(iO2) = 0.21d0

               CALL CKXTY (Xl, IWRK(ckbi), RWRK(ckbr), Yl)
               do n = 1,Nspec
                  scal(i,j,FirstSpec+n-1) = Yl(n)
               end do

               vel(i,j,1) = zero
               vel(i,j,2) =
     &              (flowrate_ox_O2+flowrate_ox_N2)/(6.d7*Pi*(R_ox_2**2-R_ox_1**2))

               if ((x.lt.x_hot).and.(y.ge.y_hot_lo).and.(y.le.y_hot_hi)) then
                  scal(i,j,Temp) = T_hot
               end if

            enddo
         enddo
      end if

      Patm = pamb / 101325.0d0
      call FORT_RHOfromPTY(lo,hi,
     &        scal(ARG_L1(scal),ARG_L2(scal),Density),  DIMS(scal),
     &        scal(ARG_L1(scal),ARG_L2(scal),Temp),     DIMS(scal),
     &        scal(ARG_L1(scal),ARG_L2(scal),FirstSpec),DIMS(scal),
     &        Patm)
      call FORT_HMIXfromTY(lo,hi,
     &        scal(ARG_L1(scal),ARG_L2(scal),RhoH),     DIMS(scal),
     &        scal(ARG_L1(scal),ARG_L2(scal),Temp),     DIMS(scal),
     &        scal(ARG_L1(scal),ARG_L2(scal),FirstSpec),DIMS(scal))

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
      integer i, j, n, zone, Y1, YN, getZone, ys, ye
      REAL_T x, y

      
      Y1 = FirstSpec + SDIM - 1
      YN = LastSpec  + SDIM - 1
      
c     Do only species equations
      if ((id .LE. YN) .AND. (id+ncomp-1.GE.Y1)) then

         ys = Max(id,Y1)
         ye = Min(id+ncomp-1,YN)
         
         if (lo(2) .EQ. domlo(2)) then
            j = lo(2)
            y = float(j)*dx(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*dx(1)+domnlo(1)
               zone = getZone(x,y)
               if (zone .EQ. 3) then
                  do n = ys - id + 1, ye - id + 1
                     diff(i,j,n) = zero
                  end do
               end if
            end do
         end if
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

      stop 'DENERROR: should no be here for H-flame problem'
      
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

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            tag(i,j) = cvmgt(set,tag(i,j),ftrac(i,j,1).gt.flametracval)
	 enddo
      enddo

      if ((level .lt. max_nozzle_lev).and.(refine_nozzle .eq. 1)) then
         do j = lo(2), hi(2)
            y = (float(j)+.5)*dx(2)+problo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*dx(1)+problo(1)
               in_refine_zone = (y - problo(2)) .le. nozzlerefheight
               in_refine_zone = in_refine_zone .and.
     &              (ABS(x-nozzlerefloc) .le. half*nozzlerefwidth)
               tag(i,j) = cvmgt(set,tag(i,j),in_refine_zone)
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
      
      return
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

#if 0
      ng = min(ARG_H1(temp)-hi(1),ARG_H2(temp)-hi(2),lo(1)-ARG_L1(temp),lo(2)-ARG_L2(temp))

      if (ng .lt. 1) then
         write(6,*) "TEMPERR cannot compute gradient, ng = ",ng
         stop
      endif

c     ::::: refine where there is temperature gradient
      if (level .lt. max_temp_lev) then
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ax = abs(temperature(i+1,j,1) - temperature(i-1,j,1))
            ay = abs(temperature(i,j+1,1) - temperature(i,j-1,1))
            aerr = max(ax,ay)
            tag(i,j) = cvmgt(set,tag(i,j),aerr.ge.tempgrad)
         enddo
      enddo
      endif
#else
      if (level .lt. max_temp_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j) = cvmgt(set,tag(i,j),temperature(i,j,1).ge.temperr)
            end do
         end do
      end if
#endif
      
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
      REAL_T    max_vort, tagval

#include "probdata.H"

      if (level .lt. max_vort_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j) = cvmgt(set,tag(i,j),ABS(vort(i,j,1)*dx(1)).ge.vorterr)
            end do
         end do
      end if

      return
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
c     rho.Yl, we need all Yls...also need to evaluate EOS since we
c     really are specifying T and Yls.  so, all this is centralized
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

      return
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

      return
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
c     rho.Yl, we need all Yls...also need to evaluate EOS since we
c     really are specifying T and Yls.  so, all this is centralized
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

      return
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
c     rho.Yl, we need all Yls...also need to evaluate EOS since we
c     really are specifying T and Yls.  so, all this is centralized
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
c     rho.Yl, we need all Yls...also need to evaluate EOS since we
c     really are specifying T and Yls.  so, all this is centralized
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

#include "cdwrk.H"
#include "bc.H"
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
c     rho.Yl, we need all Yls...also need to evaluate EOS since we
c     really are specifying T and Yls.  so, all this is centralized
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
c     rho.Yl, we need all Yls...also need to evaluate EOS since we
c     really are specifying T and Yls.  so, all this is centralized
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

      return
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
      
      subroutine FORT_ACTIVECONTROL(coft,time,dt,myproc,step,restart)

      implicit none

      REAL_T coft,time,dt,vslope,slocal,V_new,dVmax
      integer myproc,step,ierr,restart

      call bl_abort('Should not be in FORT_ACTIVECONTROL')

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
