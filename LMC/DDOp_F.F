#include "ArrayLim.H"

      subroutine DDFLUX_helper(XM,XP,YM,YP,TM,TP,CP,delta,for_T0_H1,fillAlpha,f,alpha)
      implicit none
#include "cdwrk.H"
c#include "htdata.H"
      REAL_T XM(*),XP(*),YM(*),YP(*),TM,TP,CP(*),delta,f(*),alpha(*)
      integer for_T0_H1, fillAlpha, cnt, n, m
      REAL_T Ye(maxspec), WWe, Xe(maxspec), Te, cpbe
      REAL_T dX(maxspec), dT, He(maxspec), dxInv, dx2Inv, cpf
      REAL_T rhoD(maxspec,maxspec), rhoDt(maxspec*maxspec), rhoTD(maxspec), PTC, Dmix(maxspec+1)

c     Need RU, mwt 
      REAL_T P1atm, RU, RUC, mwtInv(maxspec), pcgs, L1, L2
      integer mcdd, need_mixavg

c     Variables for evaluating coeffs at cc and moving to edges
      REAL_T CPM(maxspec), CPP(maxspec), WWM, WWP, rDmix(maxspec+1)
      REAL_T rhoDM(maxspec,maxspec), rhoDtM(maxspec*maxspec), rhoTDM(maxspec), PTCM, rDmixM(maxspec+1)
      REAL_T rhoDP(maxspec,maxspec), rhoDtP(maxspec*maxspec), rhoTDP(maxspec), PTCP, rDmixP(maxspec+1)


c     If mcdd=1, do full transport, otherwise do mixture-averaged
      mcdd = 1
      mcdd = 0
      if (fillAlpha.eq.1  .or.  mcdd.eq.0 ) then      
         need_mixavg = 1
      else
         need_mixavg = 0
      endif


      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RU = RU * 1.d-4
      call CKWT(IWRK(ckbi), RWRK(ckbr), mwtInv)
      do n=1,Nspec
         mwtInv(n) = 1.d0 / mwtInv(n)
      enddo


      dxInv = 1.d0/delta
      cpbe = 0.d0
      do n = 1,Nspec
         Ye(n)  = (  YP(n) + YM(n)  ) * 0.5d0
         Xe(n)  = (  XP(n) + XM(n)  ) * 0.5d0
         dX(n)  = (  XP(n) - XM(n)  ) * dxInv
         cpbe = cpbe + CP(n)*Ye(n)
      enddo
      Te = ( TP + TM ) * 0.5d0
      dT = ( TP - TM ) * dxInv

      CALL CKHMS(Te,IWRK(ckbi),RWRK(ckbr),He)
      CALL CKMMWY(Ye, IWRK(ckbr), RWRK(ckbr), WWe)

#if 1
      CALL CKMMWY(YM, IWRK(ckbr), RWRK(ckbr), WWM)
      CALL CKMMWY(YP, IWRK(ckbr), RWRK(ckbr), WWP)

c     Leave in cgs for eglib
      CALL CKCPMS(TM,IWRK(ckbi),RWRK(ckbr),CPM)
      CALL CKCPMS(TP,IWRK(ckbi),RWRK(ckbr),CPP)

      CALL EGSPAR(TM,XM,YM,CPM,RWRK(egbr),IWRK(egbi))
      if (mcdd.eq.1) then
         CALL EGSLTDR5(TM,YM,WWM,RWRK(egbr),IWRK(egbi),PTCM,rhoTDM,rhoDtM)
      endif
      if (need_mixavg.eq.1) then
         CALL EGSVR1(TM,YM,RWRK(egbr),rDmixM)
         CALL EGSL1( 1.d0,TM,XM,RWRK(egbr),L1)
         CALL EGSL1(-1.d0,TM,XM,RWRK(egbr),L2)
         PTCM = (L1 + L2) * .5d0
      endif

      CALL EGSPAR(TP,XP,YP,CPP,RWRK(egbr),IWRK(egbi))
      if (mcdd.eq.1) then
         CALL EGSLTDR5(TP,YP,WWP,RWRK(egbr),IWRK(egbi),PTCP,rhoTDP,rhoDtP)
      endif
      if (need_mixavg.eq.1) then
         CALL EGSVR1(TP,YP,RWRK(egbr),rDmixP)
         CALL EGSL1( 1.d0,TP,XP,RWRK(egbr),L1)
         CALL EGSL1(-1.d0,TP,XP,RWRK(egbr),L2)
         PTCP = (L1 + L2) * .5d0
      endif

      PTC = (PTCM+PTCP)*.5d-5

      if (mcdd .eq.1) then
         cnt = 1
         do n=1,Nspec
            if (Ye(n).lt.1.d-30) then
               rhoTD(n) = (rhoTDM(n)+rhoTDP(n))*.5d-1
            else
               rhoTD(n) = (YM(n)*rhoTDM(n)+YP(n)*rhoTDP(n))/Ye(n)*.5d-1
            endif
            do m=1,Nspec
               rhoD(m,n) = (rhoDtM(cnt)+rhoDtP(cnt))*.5d-1
               cnt = cnt+1
            enddo
            He(n) = He(n) * 1.d-4
         enddo
      endif

      if (need_mixavg.eq.1) then
         do n=1,Nspec
            rDmix(n) = WWe * mwtInv(n) * (rDmixM(n) + rDmixP(n)) * .5d-1
         end do
         rDmix(Nspec+1) = PTC / cpbe
      endif
#else
      CALL EGSPAR(Te,Xe,Ye,CP,RWRK(egbr),IWRK(egbi))
      if (mcdd.eq.1) then
         CALL EGSLTDR5(Te,Ye,WWe,RWRK(egbr),IWRK(egbi),PTC,rhoTD,rhoDt)
      endif
      if (need_mixavg.eq.1) then
         CALL EGSVR1(Te,Ye,RWRK(egbr),rDmix)
         CALL EGSL1( 1.d0,Te,Xe,RWRK(egbr),L1)
         CALL EGSL1(-1.d0,Te,Xe,RWRK(egbr),L2)
         PTC = (L1 + L2) * .5d0
      endif

      PTC = PTC * 1.d-5

      if (mcdd .eq.1) then
         cnt = 1
         do n=1,Nspec
            rhoTD(n) = rhoTD(n)*1.d-1
            do m=1,Nspec
               rhoD(m,n) = rhoDtM(cnt)*1.d-1
               cnt = cnt+1
            enddo
            He(n) = He(n) * 1.d-4
         enddo
      endif

      if (need_mixavg.eq.1) then
         do n=1,Nspec
            rDmix(n) = WWe * mwtInv(n) * rDmix(n) * 1.d-1
         end do
         rDmix(Nspec+1) = PTC / cpbe
      endif
#endif

      if (mcdd .eq.1) then
         f(Nspec+1) = 0.d0
         do n = 1,Nspec
            f(n) = - Ye(n) * rhoTD(n) * dT / Te
            do m = 1,Nspec
               f(n) = f(n)  -  rhoD(n,m) * dX(m)
            enddo
            f(Nspec+1) = f(Nspec+1) - (RU*Te/WWe)*rhoTD(n)*dX(n)
            if (for_T0_H1.eq.1) then
               f(Nspec+1) = f(Nspec+1) + He(n)*f(n)
            endif
         enddo
         f(Nspec+1) = f(Nspec+1)  -  PTC * dT
      else
         f(Nspec+1) = 0.d0
         do n = 1,Nspec
            f(n) = - rDmix(n) * ( YP(n) - YM(n) ) * dxInv
            if (for_T0_H1.eq.1) then
               f(Nspec+1) = f(Nspec+1) + He(n)*f(n)
            endif
         enddo
         f(Nspec+1) = f(Nspec+1)  -  PTC * dT
      endif

      if (fillAlpha.eq.1) then
         dx2Inv = dxInv*dxInv
         do n=1,Nspec+1
            alpha(n) = rDmix(n) * dx2Inv
         enddo
      endif

      if (for_T0_H1.eq.0) then
         do n=1,Nspec
            CP(n) = CP(n)*f(n)*dT
         enddo
      endif

      end
