      subroutine calc_divu(scal,divu,dx,time)
      implicit none
      include 'spec.h'
      real*8 scal(maxscal,0:nx+1)
      real*8 divu(0:nx+1)
      real*8 dx, time

      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1)
      real*8 rhoDijec(maxspec,maxspec,1:nx+1)
      real*8 rhoDiec(maxspec,1:nx+1)
      real*8 cpicc(1:maxspec,0:nx+1)
      real*8 LofS(1:maxspec+1,1:nx)
      real*8 Y(1:maxspec), mwmix, hi(1:maxspec), cpmix
      real*8 C(1:maxspec), omega(1:maxspec), rho, T
      integer n,i

      call ecCoef_and_dt(scal,PTCec,rhoTDec,rhoDijec,rhoDiec,cpicc,-1.d0,dx)
      call neg_divF_T(LofS,scal,PTCec,rhoTDec,rhoDijec,cpicc,dx)

      do i = 1,nx
         rho = scal(Density,i)
         do n=1,Nspec
            Y(n) = scal(FirstSpec+n-1,i) / rho
         enddo
         T = scal(Temp,i)
         call CKMMWY(Y,IWRK,RWRK,mwmix)
         call CKCPBS(T,Y,IWRK,RWRK,cpmix)
         call CKHMS(T,IWRK,RWRK,hi)

         call CKCPBS(T,Y,IWRK,RWRK,cpmix)
         call CKYTCP(Pcgs,T,Y,IWRK,RWRK,C)
         call CKWC(T,C,IWRK,RWRK,omega)

         divu(i) = LofS(Nspec+1,i) / (T*cpmix)
         do n=1,Nspec
            divu(i) = divu(i) + mwmix*invmwt(n)*LofS(n,i)
     &           - omega(n)*(mwt(n)*hi(n)/(cpmix*T) - mwmix)
         enddo
         divu(i) = divu(i) / rho
      enddo

c     Set grow cells
      divu(0) = divu(1)
      divu(nx+1) = divu(nx)
      end

