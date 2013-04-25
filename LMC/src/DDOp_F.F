#include <REAL.H>
#include <DDOp_F.H>

      integer function FORT_DDNCOEFS(model_DD0_MA1)
      implicit none
#include <cdwrk.H>
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
#include <cdwrk.H>
      REAL_T T, Y(*), CPMS(*), MWinv(*), rhoDi(*)
      REAL_T PTC, rhoTD(*),rhoDij(*)
      REAL_T X(maxspec), WW, Yt(maxspec), Tt, sum
      integer Full0_Mix1,n

c     HACK: Clean up state to something "reasonable"
      Tt = MAX(298.d0, MIN(6500.d0, T))
      sum = 0.d0
      do n=1,Nspec-1
         Yt(n) =  MAX(0.d0, MIN(1.d0, Y(n)))
         sum = sum + Yt(n)
      enddo
      Yt(Nspec) = 1.d0 - sum

      call CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
      call CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),WW)
      call EGSPAR(Tt,X,Yt,CPMS,EGRWRK(1),EGIWRK(1))
      CALL EGSVR1(Tt,Yt,EGRWRK(1),rhoDi)

      if (Full0_Mix1.eq.0) then
         call EGSLTDR5(Tt,Yt,WW,EGRWRK(1),EGIWRK(1),PTC,rhoTD,rhoDij)
         do n=1,Nspec
            rhoTD(n) = rhoTD(n) * 1.d-1
         end do
         do n=1,Nspec*Nspec
            rhoDij(n) = rhoDij(n) * 1.d-1
         enddo
      else
c     FIXME: Need to replace with cheaper call...leave here to debug though
         call EGSLTDR5(Tt,Yt,WW,EGRWRK(1),EGIWRK(1),PTC,rhoTD,rhoDij)
      endif
      do n=1,Nspec
         rhoDi(n) = WW * rhoDi(n) * MWinv(n) * 1.d-1
      end do
      PTC = PTC * 1.d-5
      end



