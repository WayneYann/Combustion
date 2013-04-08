
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "HEATTRANSFER_F.H"
#include "ArrayLim.H"

      subroutine FORT_SET_SCAL_NUMB(DensityIn, TempIn, TracIn, RhoHIn,
     &     FirstSpecIn, LastSpecIn)
      implicit none
#include "cdwrk.H"
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
#include "cdwrk.H"
#include "conp.H"
#include "htdata.H"
      integer nVals,n,nVals1
      REAL_T typ_vals(nVals)
      nVals1 = nVals-BL_SPACEDIM
c     Note: typical values are defaulted to zero, and may be left that way
      if (Density.gt.nVals1 
     &     .or. Temp.gt.nVals1
     &     .or. Temp.gt.nVals1
     &     .or. RhoH.gt.nVals1
     &     .or. Trac.gt.nVals1
     &     .or. LastSpec.gt.nVals) then
         call bl_pd_abort('cannot write typical values')
      endif
      typ_vals(Density+BL_SPACEDIM) = typVal_Density
      typ_vals(Temp+BL_SPACEDIM)    = typVal_Temp
      typ_vals(RhoH+BL_SPACEDIM)    = typVal_RhoH
      typ_vals(Trac+BL_SPACEDIM)    = typVal_Trac
      do n=1,Nspec
         typ_vals(FirstSpec+n-1+BL_SPACEDIM) = typVal_Y(n)
      enddo
      end

      subroutine FORT_SETTYPICALVALS(typ_vals,nVals)
#include "cdwrk.H"
#include "conp.H"
#include "htdata.H"
      integer nVals,n,nVals1
      REAL_T typ_vals(nVals)
      nVals1 = nVals-BL_SPACEDIM
c     Note: typical values are defaulted to zero, and may be left that way
      if (Density.gt.nVals1 
     &     .or. Temp.gt.nVals1
     &     .or. Temp.gt.nVals1
     &     .or. RhoH.gt.nVals1
     &     .or. Trac.gt.nVals1
     &     .or. LastSpec.gt.nVals) then
         call bl_pd_abort('cannot write typical values')
      endif
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
      integer muIsVar, lambdaIsVar, rhoDIsVar
      REAL_T muVal, lambdaVal, rhoDVal, prandtl, schmidt
      integer unityLe
#include "visc.H"
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
#include "cdwrk.H"
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
#include "cdwrk.H"
#include "conp.H"
#include "visc.H"
      data use_constant_mu     /.true./
      data use_constant_lambda /.true./
      data use_constant_rhoD   /.true./
      data constant_mu_val     /1.85e-5/
      data constant_lambda_val /1.85e-5/
      data constant_rhoD_val   /1.85e-5/
      data Pr / one /
      data Sc / one /
      data thickFacTR / one /

      data typVal_Density / 0.d0 /
      data typVal_Temp / 0.d0 /
      data typVal_RhoH / 0.d0 /
      data typVal_Trac / 0.d0 /
      data typVal_Y / maxspec*0.d0 /

      end


      subroutine FORT_SET_COMMON(time1,iteration1)
      REAL_T time1
      integer iteration1
#include "timedata.H"

      time = time1
      iteration = iteration1
      return
      end
      
      

