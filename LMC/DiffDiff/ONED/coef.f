      subroutine ecCoef_and_dt(S,PTCec,rhoTDec,rhoDijec,rhoDiec,cpicc,dt,dx)
      implicit none
      include 'spec.h'

      real*8  S(maxscal,0:nx+1)

      real*8  PTCec(1:nx+1)
      real*8  rhoTDec(maxspec,1:nx+1)
      real*8  rhoDijec(maxspec,maxspec,1:nx+1)
      real*8  rhoDiec(maxspec,1:nx+1)
      real*8  cpicc(1:maxspec,0:nx+1)

      real*8  PTCcc(0:nx+1)
      real*8  rhoTDcc(maxspec,0:nx+1)
      real*8  rhoDijcc(maxspec,maxspec,0:nx+1)
      real*8  rhoDicc(maxspec,0:nx+1)
      real*8  LofS(maxspec+1,1:nx)

      real*8 dx, dt, mass(maxspec)
      real*8 rhom, rhop, rhoc, ym, yp, yc, cpb
      integer i, n, m

      real*8 T

c     Compute cc transport coeffs
      do i=0,nx+1
         do n=1,Nspec
            mass(n) = S(FirstSpec+n-1,i) / S(Density,i)
         enddo
         call calc_beta(S(Temp,i),mass,S(Density,i),
     &        PTCcc(i),rhoTDcc(1,i),rhoDijcc(1,1,i),rhoDicc(1,i),cpicc(1,i))
      enddo

c     If requested, compute timestep based on Di,m and lambda/(rho.cpb)
      if (dt.gt.0) then
         dt = big
         do i=0,nx+1
            do n=1,Nspec
               if (rhoDicc(n,i) .gt. small) then
                  dt = MIN(dt,dtRedFac*S(FirstSpec+n-1,i)*S(Density,i)*dx*dx/(2.d0*rhoDicc(n,i)))
               endif
            enddo
            if (PTCcc(i) .gt. small) then
               cpb = 0.d0
               do n=1,Nspec
                  cpb = cpb + cpicc(n,i) * S(FirstSpec+n-1,i) / S(Density,i)
               enddo
               dt = MIN(dt,dtRedFac*dx*dx*S(Density,i)*cpb/(2.d0*PTCcc(i)))
            endif
         enddo
         
         if (dt.le.smallDt) then
            print *,'dt too small',dt
            stop
         endif
      endif
      
c     Compute ec transport coeffs
      do i=1,nx+1
         PTCec(i) = 0.5d0 * ( PTCcc(i-1) + PTCcc(i) )
         rhop = S(Density,i)
         rhom = S(Density,i-1)
         rhoc = 0.5d0*(rhop+rhom)
         do n=1,Nspec
            yp = S(FirstSpec+n-1,i)/S(Density,i)
            ym = S(FirstSpec+n-1,i-1)/S(Density,i-1)
            yc = 0.5d0*(yp+ym)
            if (yc.lt.small) then
               rhoTDec(n,i) = 0.5d0*(rhoTDcc(n,i-1)+rhoTDcc(n,i))
            else
               rhoTDec(n,i) = 0.5d0*(ym*rhoTDcc(n,i-1)+yp*rhoTDcc(n,i))/yc
            endif
            rhoDiec(n,i) = 0.5d0*(rhoDicc(n,i-1)+rhoDicc(n,i))
            do m=1,Nspec
               rhoDijec(n,m,i) = 0.5d0*(rhoDijcc(n,m,i-1)+rhoDijcc(n,m,i))
            enddo
         enddo         
      enddo
      end

      subroutine calc_beta(T,Y,rho,PTC,rhoTD,rhoDij,rhoDi,CPMS)
      include 'spec.h'
      real*8 T, Y(maxspec), rho, CPMS(maxspec), X(maxspec), WW, PTC
      real*8 rhoTD(maxspec),rhoDij(maxspec,maxspec)
      real*8 rhoDijt(maxspec*maxspec), rhoDi(maxspec)
      integer n,m,cnt

      call CKYTX(Y,IWRK,RWRK,X)
      call CKMMWY(Y,IWRK,RWRK,WW)
      call CKCPMS(T,IWRK,RWRK,CPMS)
      call EGSPAR(T,X,Y,CPMS,EGRWRK,EGIWRK)

      if (rhoInTrans.eq.1) then
         call EGSLTDR5(T,Y,WW,EGRWRK,EGIWRK,PTC,rhoTD,rhoDijt)
         cnt = 1
         do n=1,Nspec
            do m=1,Nspec
               rhoDij(m,n) = rhoDijt(cnt)
               cnt = cnt+1
            enddo
         enddo
         CALL EGSVR1(T,Y,EGRWRK,rhoDi)
         do n=1,Nspec
            rhoDi(n) = WW * rhoDi(n) / mwt(n)
         end do
      else
         call EGSLTD5(Pcgs,T,Y,WW,EGRWRK,EGIWRK,PTC,rhoTD,rhoDijt)
         rhoTD = rho * rhoTD
         cnt = 1
         do n=1,Nspec
            do m=1,Nspec
               rhoDij(m,n) = rhoDijt(cnt)
               cnt = cnt+1
            enddo
         enddo
         CALL EGSV1(Pcgs,T,Y,WW,EGRWRK,rhoDi)
         do n=1,Nspec
            rhoDi(n) = rho * WW * rhoDi(n) / mwt(n)
         end do
      endif
      end

