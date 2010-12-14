#include "DDOp_F.H"

      integer function FORT_DDNCOEFS(model_DD0_MA1)
      implicit none
#include "cdwrk.H"
      integer model_DD0_MA1,NfullDD, NmixAvg

      DMIXoffset = 0
      PTCoffset = DMIXoffset + Nspec
      NmixAvg = PTCoffset + 1
      
      TDoffset = NmixAvg
      Doffset = TDoffset + Nspec
      NfullDD = Doffset + Nspec*Nspec

      if (model_DD0_MA1.eq.0) then
         Ncoefs = NfullDD
      else if (model_DD0_MA1.eq.1) then
         Ncoefs = NmixAvg
      else if (model_DD0_MA1.lt.0 .or. model_DD0_MA1.gt.1) then
         print *,'FORT_DDNCOEFS: invalid model'
         call bl_pd_abort()
      endif
      FORT_DDNCOEFS = Ncoefs
      end

      subroutine calc_coefs(T,Y,CPMS,MWinv,rhoDi,Full0_Mix1,PTC,rhoTD,rhoDij)
      implicit none
#include "cdwrk.H"
      REAL_T T, Y(*), CPMS(*), MWinv(*), rhoDi(*)
      REAL_T PTC, rhoTD(*),rhoDij(*)
      REAL_T X(maxspec), WW
      integer Full0_Mix1,n

      call CKYTX(Y,IWRK(ckbi),RWRK(ckbr),X)
      call CKMMWY(Y,IWRK(ckbi),RWRK(ckbr),WW)
      call EGSPAR(T,X,Y,CPMS,RWRK(egbr),IWRK(egbi))
      CALL EGSVR1(T,Y,RWRK(egbr),rhoDi)
      if (Full0_Mix1.eq.0) then
         call EGSLTDR5(T,Y,WW,RWRK(egbr),IWRK(egbi),PTC,rhoTD,rhoDij)
         do n=1,Nspec
            rhoTD(n) = rhoTD(n) * 1.d-1
         end do
         do n=1,Nspec*Nspec
            rhoDij(n) = rhoDij(n) * 1.d-1
         enddo
      else
c     FIXME: Need to replace with cheaper call...leave here to debug though
         call EGSLTDR5(T,Y,WW,RWRK(egbr),IWRK(egbi),PTC,rhoTD,rhoDij)
      endif
      do n=1,Nspec
         rhoDi(n) = WW * rhoDi(n) * MWinv(n) * 1.d-1
      end do
      PTC = PTC * 1.d-5
      end



