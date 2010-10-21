      subroutine lmc()
      implicit none
      include 'spec.h'
      integer nsteps

      real*8  scal_old(maxscal,0:nx+1)
      real*8  scal_new(maxscal,0:nx+1)
      real*8  PTCcc(0:nx+1)
      real*8  rhoTDcc(maxspec,0:nx+1)
      real*8  rhoDijcc(maxspec,maxspec,0:nx+1)
      real*8  rhoDicc(maxspec,0:nx+1)
      real*8  Y(maxspec,0:nx+1)
      real*8  LofS(maxspec+1,1:nx)

      real*8  PTCec_old(1:nx+1)
      real*8  rhoTDec_old(maxspec,1:nx+1)
      real*8  rhoDijec_old(maxspec,maxspec,1:nx+1)
      real*8  rhoDiec_old(maxspec,1:nx+1)

      real*8 dx, dt, problo, probhi, enth, cpb, big, small, smallDt
      real*8 rhom, rhop, rhoc, ym, yp, yc
      real*8 sum, maxsum, maxDiff, avgMag, T, rho, Yhalf
      integer i, Npmf, n, m
      real*8 x, time
      real*8 Patm, flame_offset, pmfdata(maxspec+3), mole(maxspec), mass(maxspec)
      real*8 dtRedFac

      integer Niter, maxIters, step
      integer probtype, alt_spec_update
      real*8 res(NiterMAX)

c     Initialize chem/tran database
      call initchem()

c     Set defaults, change with namelist
      nsteps = 1000
      problo = 0.0
      probhi = 3.5
      flame_offset = 0.d0
      Patm = 1.d0
      probtype = 2
      dtRedFac = 1.d0
      big = 1.d30
      small = 1.d-30
      smallDt = 1.d-12
      alt_spec_update = 0

      call CKRP(IWRK,RWRK,RU,RUC,P1ATM)
      Pcgs = Patm * P1ATM

      Density = 1
      FirstSpec = Density + 1
      LastSpec = FirstSpec + Nspec - 1
      RhoH = LastSpec + 1
      Temp = RhoH + 1
      dx = (probhi-problo)/nx

      do i=1,nx
         x = (i+0.5d0)*dx - flame_offset
         if (probtype.eq.1) then
            call pmf(x,x,pmfdata,Npmf)
            if (Npmf.ne.Nspec+3) then
               print *,'mismatched pmf'
               stop
            endif
            scal_new(Temp,i) = pmfdata(1)
            do n=1,Nspec
               mole(n) = pmfdata(3+n)
            enddo
         else
            scal_new(Temp,i) = 298.d0
            do n=1,Nspec
               mole(n) = 0.d0
            enddo
            if (x.lt.0.5d0*(problo+probhi)) then
               mole(1) = 0.3d0
               mole(4) = 0.3d0
               mole(9) = 0.4d0
            else
               mole(4) = 0.21d0
               mole(9) = 0.79d0
            endif
         endif

         call CKXTY(mole,IWRK,RWRK,mass)         
         call CKRHOY(Pcgs,scal_new(Temp,i),mass,IWRK,RWRK,scal_new(Density,i))
         call CKHBMS(scal_new(Temp,i),mass,IWRK,RWRK,scal_new(RhoH,i))

         do n=1,Nspec
            scal_new(FirstSpec+n-1,i) = mass(n) * scal_new(Density,i)
         enddo
         scal_new(RhoH,i) = scal_new(RhoH,i) * scal_new(Density,i)
         
      enddo
      
c     Left boundary grow cell
      scal_new(Density,0) = scal_new(Density,1)
      scal_new(Temp,0) = scal_new(Temp,1)
      do n=1,Nspec
         scal_new(FirstSpec+n-1,0) = scal_new(FirstSpec+n-1,1)
      enddo
      scal_new(RhoH,0) = scal_new(RhoH,1)

c     Right boundary grow cell
      scal_new(Density,nx+1) = scal_new(Density,nx)
      scal_new(Temp,nx+1) = scal_new(Temp,nx)
      do n=1,Nspec
         scal_new(FirstSpec+n-1,nx+1) = scal_new(FirstSpec+n-1,nx)
      enddo
      scal_new(RhoH,nx+1) = scal_new(RhoH,nx)



      time = 0.d0
      call print_soln(time,scal_new,'soln_start.dat',dx,problo)
      do step=1,nsteps

c     Update state and compute Yold
         do i=0,nx+1
            do n=1,maxscal
               scal_old(n,i) = scal_new(n,i)
            enddo
         enddo
         
         do i=0,nx+1
            do n=1,Nspec
               Y(n,i) = scal_old(FirstSpec+n-1,i) / scal_old(Density,i)
            enddo
         enddo

c     Compute cc transport coeffs
         do i=0,nx+1
            call calc_beta(scal_old(Temp,i),Y(1,i),
     &           PTCcc(i),rhoTDcc(1,i),rhoDijcc(1,1,i),rhoDicc(1,i))
         enddo

c     Compute timestep based on Di,m and lambda/(rho.cpb)
         dt = big
         do i=0,nx+1
            do n=1,Nspec
               if (rhoDicc(n,i) .gt. small) then
                  dt = MAX(small,
     &                     MIN(dt,dtRedFac*scal_old(FirstSpec+n-1,i)*dx*dx/(2.d0*rhoDicc(n,i))))
               endif
            enddo
            call CKCPBS(scal_old(Temp,i),mass,IWRK,RWRK,cpb)
            if (PTCcc(i) .gt. small) then
               dt = MAX(small,
     &                  MIN(dt,dtRedFac*dx*dx*scal_old(Density,i)*cpb/(2.d0*PTCcc(i))))
            endif

         enddo
         if (dt.lt.smallDt) then
               print *,'dt too small at i=',i,': ',dt
               goto 100
         endif

c     Compute ec transport coeffs
         do i=1,nx+1
            PTCec_old(i) = 0.5d0 * ( PTCcc(i-1) + PTCcc(i) )
            rhop = scal_old(Density,i)
            rhom = scal_old(Density,i-1)
            rhoc = 0.5d0*(rhop+rhom)
            do n=1,Nspec
               yp = scal_old(FirstSpec+n-1,i)/rhop
               ym = scal_old(FirstSpec+n-1,i-1)/rhom
               yc = 0.5d0*(yp+ym)
               if (yc.lt.small) then
                  rhoTDec_old(n,i) = 0.5d0*(rhoTDcc(n,i-1)+rhoTDcc(n,i))
               else
                  rhoTDec_old(n,i) = 0.5d0*(ym*rhoTDcc(n,i-1)/rhom+yp*rhoTDcc(n,i)/rhop)*rhoc/yc
               endif
               rhoDiec_old(n,i) = 0.5d0*(rhoDicc(n,i-1)/rhom+rhoDicc(n,i)/rhop)*rhoc
               do m=1,Nspec
                  rhoDijec_old(n,m,i) = 0.5d0*(rhoDijcc(n,m,i-1)+rhoDijcc(n,m,i))
               enddo
            enddo
         enddo

         maxsum = 0.d0
         do i=1,nx
            sum = 0.d0
            do n=1,Nspec
               sum = sum + rhoTDcc(n,i)*scal_old(FirstSpec+n-1,i)/scal_old(Density,i)**2
            enddo
            maxsum = MAX(maxsum,ABS(sum))
         enddo
         print *,'  cc max sum:',maxsum
         
         maxsum = 0.d0
         do i=1,nx+1
            rhop = scal_old(Density,i)
            rhom = scal_old(Density,i-1)
            rhoc = 0.5d0*(rhop+rhom)
            sum = 0.d0
            do n=1,Nspec
               yp = scal_old(FirstSpec+n-1,i)/rhop
               ym = scal_old(FirstSpec+n-1,i-1)/rhom
               yc = 0.5d0*(yp+ym)
c               sum = sum + yc*rhoTDec_old(n,i)/rhoc
               do m=1,Nspec
                  sum = sum + rhoDijec_old(n,m,i)
               enddo
            enddo
            maxsum = MAX(ABS(sum),maxsum)
         enddo
         print *,'  ec max sum:',maxsum
         

         
         call LinOpApply(LofS,scal_old,PTCec_old,rhoTDec_old,rhoDijec_old,dx)

c     Form explicit update
         do i=1,nx+1

            if (alt_spec_update .eq. 0) then

               do n=1,Nspec
                  scal_new(FirstSpec+n-1,i)=scal_old(FirstSpec+n-1,i) +
     &                 dt*LofS(n,i)
               enddo
               scal_new(Density,i) = scal_old(Density,i)

            else if (alt_spec_update .eq. 1) then

               sum = 0.d0
               do n=1,Nspec-1
                  scal_new(FirstSpec+n-1,i)=scal_old(FirstSpec+n-1,i) +
     &                 dt*LofS(n,i)
                  sum = sum + scal_new(FirstSpec+n-1,i)
               enddo
               scal_new(Density,i) = scal_old(Density,i)
               scal_new(LastSpec,i) = scal_new(Density,i) - sum

            else if (alt_spec_update .eq. 2) then

               sum = 0.d0
               do n=1,Nspec
                  scal_new(FirstSpec+n-1,i)=scal_old(FirstSpec+n-1,i) +
     &                 dt*LofS(n,i)
                  sum = sum + scal_new(FirstSpec+n-1,i)
               enddo
               scal_new(Density,i) = sum

            else
               print *,'invalid value for alt_spec_update: ',alt_spec_update
            endif

            scal_new(RhoH,i)=scal_old(RhoH,i) + dt*LofS(Nspec+1,i)
         enddo
         
c     Recompute temperature
         maxIters=0
         do i=1,nx
            do n=1,Nspec
               mass(n) = scal_new(FirstSpec+n-1,i)/scal_new(Density,i)
            enddo
            enth = scal_new(RhoH,i)/scal_new(Density,i)
            call FORT_TfromHYpt(scal_new(Temp,i),enth,mass,
     &           Nspec,errMax,NiterMAX,res,Niter)
            if (Niter.lt.0) then
               print *,'RhoH->T failed at i=',i
               goto 100
            endif
            maxIters = MAX(maxIters,Niter)
         enddo

c         print *,'t=',time,' RhoH->T maxIters: ',maxIters
         time = time + dt
         print *,'step=', step, ' t=',time,' dt=',dt
      enddo
 100  call print_soln(time,scal_new,'soln_end.dat',dx,problo)

      end

      subroutine LinOpApply(LofS,S,PTCec,rhoTDec,rhoDijec,dx)
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1),rhoDijec(maxspec,maxspec,1:nx+1),dx
      real*8 coef(maxspec+1,1:nx)
      real*8 Y(maxspec,0:nx+1), X(maxspec,0:nx+1), WWe, CPMS,PTC
      real*8 de(maxspec+1,1:nx+1), q(1:nx+1), F(maxspec,1:nx+1)
      real*8 Ye(maxspec), Te, dxInv2, He(maxspec)
      real*8 sum, maxsum, rhoe, enthe
      integer i,n,m,setTfromH,Niter,maxIter
      real*8 res(NiterMAX)

      setTfromH = 0
      dxInv2 = 1.d0/(dx*dx)
      maxsum=0.d0
      do i=0,nx+1
         sum = 0.d0
         do n=1,Nspec
            Y(n,i) = S(FirstSpec+n-1,i) / S(Density,i) 
            sum = sum + Y(n,i)
         enddo
         maxsum = MAX(ABS(1.d0 - sum),maxsum)
         call CKYTX(Y(1,i),IWRK,RWRK,X(1,i))
      enddo
      print *,'LinOp: maxsum Y = ',maxsum

      maxsum=0.d0
      maxiter=0
      do i=1,nx+1
         do n=1,Nspec
            de(n,i) = (X(n,i)-X(n,i-1)) / dx
         enddo
         de(Nspec+1,i) = (S(Temp,i)-S(Temp,i-1)) / dx

         do n=1,Nspec
            Ye(n) = 0.5d0*(Y(n,i)+Y(n,i-1))
         enddo
         call CKMMWY(Ye,IWRK,RWRK,WWe)

         Te = 0.5d0*(S(Temp,i)+S(Temp,i-1))
         rhoe = 0.5d0*(S(Density,i)+S(Density,i-1))
         if (setTfromH.eq.1) then
            enthe = 0.5d0*(S(RhoH,i)+S(RhoH,i-1))/rhoe
            call FORT_TfromHYpt(Te,enthe,Ye,Nspec,errMax,NiterMAX,res,Niter)
            if (Niter.lt.0) then
               print *,'RhoH->T failed at i=',i
               stop
            endif
            maxiter=MAX(maxiter,Niter)
         else if (setTfromH.eq.0) then
            Te = Pcgs * WWe / (rhoe * RU)
         endif
         call CKHMS(Te,IWRK,RWRK,He) 

         q(i) = 0.d0
         sum = 0.d0
         do n = 1,Nspec
            F(n,i) = - Ye(n)*rhoTDec(n,i)*de(Nspec+1,i)/Te
            do m = 1,Nspec
               F(n,i) = F(n,i) - rhoDijec(n,m,i)*de(m,i)
            enddo
c            sum = sum + F(n,i)
            q(i) = q(i) - (RU*Te/WWe)*rhoTDec(n,i)*de(n,i) + He(n)*F(n,i)
         enddo
         sum = (Pcgs - rhoe*RU*Te/WWe)/Pcgs
         maxsum = MAX(ABS(sum),maxsum)
         q(i) = q(i) - PTCec(i)*de(Nspec+1,i)
      enddo
      print *,'LinOp: maxsum Y.theta = ',maxsum
      if (setTfromH.eq.1) then
         print *,'    maxiter: ',maxiter
      endif

      maxsum=0.d0
      do i=1,nx
         sum = 0.d0
         do n = 1,Nspec
            LofS(n,i) = - dxInv2*(F(n,i+1) - F(n,i))
            sum = sum + LofS(n,i)
         enddo
         maxsum = MAX(sum,maxsum)
         LofS(Nspec+1,i) = - dxInv2*(q(i+1) - q(i))
      enddo
      print *,'LinOp: maxsum L = ',maxsum*dx*dx

      end


      subroutine calc_beta(T,Y,PTC,rhoTD,rhoDij,rhoDi)
      include 'spec.h'
      real*8 T, Y(maxspec), X(maxspec), WW, CPMS,PTC
      real*8 rhoTD(maxspec),rhoDij(maxspec,maxspec),rhoDijt(maxspec*maxspec)
      real*8 rhoDi(maxspec)
      integer n,m,cnt

      call CKYTX(Y,IWRK,RWRK,X)
      call CKMMWY(Y,IWRK,RWRK,WW)
      call CKCPBS(T,Y,IWRK,RWRK,CPMS)
      call EGSPAR(T,X,Y,CPMS,EGRWRK,EGIWRK)
      call EGSLTDR5(T,Y,WW,EGRWRK,EGIWRK,PTC,rhoTD,rhoDijt)
      cnt = 1
      do n=1,Nspec
         do m=1,Nspec
            rhoDij(m,n) = rhoDijt(cnt)
            cnt = cnt+1
         enddo
      enddo

c     Mixture-averaged transport coefficients
      CALL EGSV1(Pcgs,T,Y,WW,EGRWRK,rhoDi)
      do n=1,Nspec
         rhoDi(n) = Y(n) * WW * rhoDi(n) / mwt(n)
      end do
      end

      integer function FORT_GETCKSPECNAME(i, coded)
      include 'spec.h'
      integer i
      integer coded(*)
      integer names(maxspec*maxspnml)
      integer j, str_len
      call CKSYMS(names, maxspnml)
      do j = 1, maxspnml
         coded(j) = names(maxspnml*(i-1)+j)
      end do
      do j = 1, maxspnml
         if (coded(j).eq.ICHAR(' ')) then
            str_len = j
            exit
         endif 
      end do
      FORT_GETCKSPECNAME = str_len - 1
      end

      integer function get_spec_name(name, j)
      include 'spec.h'
      integer i, j, FORT_GETCKSPECNAME
      integer coded(maxspnml)
      character*(maxspnml) name
      get_spec_name = FORT_GETCKSPECNAME(j, coded)
      do i = 1, maxspnml
         name(i:i) = ' '
      end do
      do i = 1, get_spec_name
         name(i:i) = char(coded(i))
      end do
      end


      subroutine print_soln(time,scal,filename,dx,plo)
      include 'spec.h'
      real*8 time, scal(maxscal,0:nx+1), dx, plo
      character*(*) filename
      character*(maxspnml) names(maxspec)
      integer n,i,j, get_spec_name, nlen(maxspec)
      do n=1,Nspec
         nlen = get_spec_name(names(n), n)
      enddo
      open(unit=12,file=filename)
      write(12,'(50a)') 'VARIABLES=X Rho ',(names(n),n=1,Nspec),
     &     ' RhoH Temp'
      write(12,'(a,i5,a,g20.8,a)') 'ZONE I=',nx,' T= "',time,
     &     '" DATAPACKING=POINT'
      do i=1,nx
         write(12,'(50g20.8)') (i+0.5d0)*dx + plo, scal(1,i),
     &        (scal(1+n,i)/scal(1,i),n=1,Nspec),scal(Nspec+2,i),scal(Nspec+3,i)
      enddo
      close(12)
      end

