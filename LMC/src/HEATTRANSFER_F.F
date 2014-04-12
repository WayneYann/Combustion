
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include <BC_TYPES.H>
#include <HEATTRANSFER_F.H>
#include <ArrayLim.H>

      subroutine FORT_SET_SCAL_NUMB(DensityIn, TempIn, TracIn, RhoHIn,
     &     FirstSpecIn, LastSpecIn)
      implicit none
#include <cdwrk.H>
#include <htdata.H>
      integer DensityIn, TempIn, TracIn, RhoHIn, FirstSpecIn, LastSpecIn
c
c ::: Remove SPACEDIM from the counter, since those spots contain the
c ::: velocity, and our INITDATA function below fills the scalar state
c ::: However, add one since the C++ is 0-based      
c     
      Density = DensityIn - BL_SPACEDIM + 1
      Temp = TempIn - BL_SPACEDIM + 1
      Trac = TracIn - BL_SPACEDIM + 1
      RhoH = RhoHIn - BL_SPACEDIM + 1
      FirstSpec = FirstSpecIn - BL_SPACEDIM + 1
      LastSpec = LastSpecIn - BL_SPACEDIM + 1
      end      
      
      subroutine FORT_GETTYPICALVALS(typ_vals,nVals)
      implicit none
#include <cdwrk.H>
#include <conp.H>
#include <htdata.H>
      integer nVals,n,nVals1
      REAL_T typ_vals(nVals)
      nVals1 = nVals-BL_SPACEDIM
c     Note: typical values are defaulted to zero, and may be left that way
      if (Density.gt.nVals1 
     &     .or. Temp.gt.nVals1
     &     .or. RhoH.gt.nVals1
     &     .or. Trac.gt.nVals1
     &     .or. LastSpec.gt.nVals) then
         call bl_pd_abort('cannot write typical values')
      endif
      do n=1,BL_SPACEDIM
         typ_vals(n) = typVal_Vel
      enddo
      typ_vals(Density+BL_SPACEDIM) = typVal_Density
      typ_vals(Temp+BL_SPACEDIM)    = typVal_Temp
      typ_vals(RhoH+BL_SPACEDIM)    = typVal_RhoH
      typ_vals(Trac+BL_SPACEDIM)    = typVal_Trac
      do n=1,Nspec
         typ_vals(FirstSpec+n-1+BL_SPACEDIM) = typVal_Y(n)
      enddo
      end


      subroutine FORT_SETTYPICALVALS(typ_vals,nVals)
      implicit none
#include <cdwrk.H>
#include <conp.H>
#include <htdata.H>
      integer nVals,n,nVals1
      REAL_T typ_vals(nVals)
      nVals1 = nVals-BL_SPACEDIM
c     Note: typical values are defaulted to zero, and may be left that way
      if (Density.gt.nVals1 
     &     .or. Temp.gt.nVals1
     &     .or. RhoH.gt.nVals1
     &     .or. Trac.gt.nVals1
     &     .or. LastSpec.gt.nVals) then
         call bl_pd_abort('cannot write typical values')
      endif
      do n=1,BL_SPACEDIM
         typVal_Vel = typ_vals(n)
      enddo
      typVal_Density = typ_vals(Density+BL_SPACEDIM)
      typVal_Temp    = typ_vals(Temp+BL_SPACEDIM)
      typVal_RhoH    = typ_vals(RhoH+BL_SPACEDIM)
      typVal_Trac    = typ_vals(Trac+BL_SPACEDIM)
      do n=1,Nspec
         typVal_Y(n) = typ_vals(FirstSpec+n-1+BL_SPACEDIM)
      enddo
      end


      subroutine FORT_SET_HT_VISC_COMMON(muIsVar,     muVal,
     &                                   lambdaIsVar, lambdaVal,
     &                                   rhoDIsVar,   rhoDVal,
     &                                   prandtl, schmidt, unityLe)
      implicit none
      integer muIsVar, lambdaIsVar, rhoDIsVar
      REAL_T muVal, lambdaVal, rhoDVal, prandtl, schmidt
      integer unityLe
#include <visc.H>
      if (muIsVar .EQ. 1) then
         use_constant_mu = .false.
         constant_mu_val = - one
      else
         use_constant_mu = .true.
         constant_mu_val = muVal
      end if
      if (lambdaIsVar .EQ. 1) then
         use_constant_lambda = .false.
         constant_lambda_val = - one
      else
         use_constant_lambda = .true.
         constant_lambda_val = lambdaVal
      end if
      if (rhoDIsVar .EQ. 1) then
         use_constant_rhoD = .false.
         constant_rhoD_val = - one
      else
         use_constant_rhoD = .true.
         constant_rhoD_val = rhoDVal
      end if
      Pr = prandtl
      Sc = schmidt
      LeEQ1 = unityLe .ne. 0
      end
      
c-----------------------------------------------------------------------

      subroutine FORT_GETPAMB(pambout, dpdt_factorout)
      implicit none
#include <cdwrk.H>
#include <htdata.H>
      REAL_T dpdt_factorout, pambout

      dpdt_factorout = dpdt_factor
      pambout = pamb
      return
      end

c=======================================================================

      block data FORT_INIT_VISC_DATA

c ... what the hell kind of dinosuar is this? Well, block
c     data is the only guranteed way (i.e., portable)
c     to use data statements to initialize common block variables.
c     block data is almost certainly not allowed in f90, but this
c     routine should not be too painful to reimplement. I would just
c     suck its functionality into some other initialization routine.
c     -- rbp

      implicit none
#include <cdwrk.H>
#include <conp.H>
#include <visc.H>
      data use_constant_mu     /.true./
      data use_constant_lambda /.true./
      data use_constant_rhoD   /.true./
      data constant_mu_val     /1.85d-5/
      data constant_lambda_val /1.85d-5/
      data constant_rhoD_val   /1.85d-5/
      data Pr / one /
      data Sc / one /
      data thickFacTR / one /

      data typVal_Density / 0.d0 /
      data typVal_Temp / 0.d0 /
      data typVal_RhoH / 0.d0 /
      data typVal_Trac / 0.d0 /
      data typVal_Y / maxspec*0.d0 /
      data typVal_YMAX / 1.d0 /
      data typVal_YMIN / 1.d-6 /

      end


      subroutine FORT_SET_COMMON(time1,iteration1)
      implicit none
      REAL_T time1
      integer iteration1
#include <timedata.H>

      time = time1
      iteration = iteration1
      return
      end

      subroutine FORT_ACTIVECONTROL(coft,time,dt,myproc,step,restart,usetemp)

      implicit none

#include <probdata.H>
#include <cdwrk.H>
#include <bc.H>

      REAL_T coft,time,dt,vslope,slocal,V_new,dVmax,dVmin
      integer myproc,step,ierr,restart,usetemp
      REAL_T r1,r2,r3,r4,r5,r6,r7
      REAL_T alpha,xsmb,vpmax,exp1 
      REAL_T rhs1,rhs2,vt_tay,vtt_tay
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

      if (usetemp.eq.0 .and. coft_old .lt. zero) coft_old = coft

      zbase_control = zbase_control + V_in*dt + dV_control*dt**2
      V_in_old = V_in
      V_in = V_in + dt*dV_control

      slocal = half*(V_in_old + V_in) - (coft - coft_old)/(dt*scale_control)
      sest = (one - corr)*sest + corr*slocal
#if 1
c    linear
c     vslope = two*((cfix-coft)/(scale_control*tau_control) + sest - V_in)/tau_control

c     V_new = V_in + dt*vslope

c    quadratic 1
c     vslope = 3.d0*((cfix-coft)/(scale_control*tau_control) + sest - V_in)/tau_control
c     V_new = V_in + (dt-0.5d0*dt**2/tau_control)*vslope

c    quadratic 2
      rhs2 = 2.d0*((cfix-coft)/(scale_control*tau_control) + sest - V_in)/tau_control
      rhs1 = (sest - V_in)/tau_control

      vt_tay = 3.d0*rhs2 - 2.d0*rhs1
      vtt_tay = 6.d0*(rhs1-rhs2)/tau_control
      V_new = V_in + dt*vt_tay + 0.5d0*dt**2*vtt_tay


      dVmax = changeMax_control * one
      dVmin = changeMax_control * max(one,V_in)
#else


      vpmax = max(V_in,V_in_old,one)*changeMax_control/dt
      xsmb = -(cfix-coft)/scale_control

      alpha = sqrt(abs(xsmb)/vpmax)

      exp1 = dexp(-2.d0*tau_control/alpha)

      V_new = sest+ (xsmb)/alpha*
     &   (exp1*dexp(dt/alpha)-dexp(-dt/alpha))/
     &      (1.d0+exp1)


      dVmax = changeMax_control * max(one,V_in)
#endif
      V_new = MIN(MAX(V_new,V_in-dVmin),V_in+dVmax)
      V_new = MAX(zero,V_new)
      V_new = Min(V_new,controlVelMax)

      tbase_control = time
      dV_control = (V_new - V_in)/dt

      if (myproc.eq.0) then
         print *
         print *,'****************** control:', scale_control
         print *,'time,dt,coft,V_new:',time,dt,coft,V_new
         print *,'changeMax_control: ', changeMax_control
#if 1
         print *,'vslope,sest,cfix,slocal:',vslope,sest,cfix,slocal
#else
         print *,'alpha,sest,cfix,slocal:',alpha,sest,cfix,slocal
#endif
         print *,'coft_old,V_in_old:',coft_old,V_in_old
         call flush(6)
      endif

      coft_old = coft

      if (myproc.eq.0) then
         open(13,file=ac_hist_file,form='formatted',position='append')
         write(13,1000) step,time,V_in,tbase_control,zbase_control,
     &        dV_control,sest,coft_old
         close(13)
      endif
 1000 format(i7,7g26.18)
      end
