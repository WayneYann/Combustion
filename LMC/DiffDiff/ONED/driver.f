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

      real*8  LofS(maxspec+1,1:nx)
      real*8  PTCec_new(1:nx+1)
      real*8  rhoTDec_new(maxspec,1:nx+1)
      real*8  rhoDijec_new(maxspec,maxspec,1:nx+1)
      real*8  rhoDiec_new(maxspec,1:nx+1)

      real*8  PTCec_old(1:nx+1)
      real*8  rhoTDec_old(maxspec,1:nx+1)
      real*8  rhoDijec_old(maxspec,maxspec,1:nx+1)
      real*8  rhoDiec_old(maxspec,1:nx+1)

      real*8 dx, dt, problo, probhi, enth
      integer i, k, Npmf, n, m
      real*8 x, time
      real*8 Patm, flame_offset, pmfdata(maxspec+3), mole(maxspec), mass(maxspec)
      real*8 redFac

      integer Niter, NiterMAX, maxIters
      parameter (NiterMAX=20)
      real*8 res(NiterMAX), errMAX

c     Initialize chem/tran database
      call initchem()

c     Set defaults, change with namelist
      nsteps = 1
      problo = 0.0
      probhi = 3.5
      flame_offset = 0.d0
      Patm = 1.d0
      errMAX = 1.d-8

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
c         call pmf(x,x,pmfdata,Npmf)
c         if (Npmf.ne.Nspec+3) then
c            print *,'mismatched pmf'
c            stop
c         endif
c         scal_new(Temp,i) = pmfdata(1)
c         do n=1,Nspec
c            mole(n) = pmfdata(3+n)
c         enddo
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
      do k=1,nsteps
         do i=0,nx+1
            do n=1,maxscal
               scal_old(n,i) = scal_new(n,i)
            enddo
         enddo
         
c     Compute cc transport coeffs
         redFac = 1.d-4
         dt = 1.e30
         do i=0,nx+1
            call calc_beta(scal_old(Density,i),scal_old(Temp,i),scal_old(FirstSpec,i),
     &           PTCcc(i),rhoTDcc(1,i),rhoDijcc(1,1,i),rhoDicc(1,i))
            do n=1,Nspec
               dt = MIN(dt,redFac*0.5d0*rhoDicc(n,i)/(dx*dx))
            enddo
         enddo
         
c     Compute ec transport coeffs
         do i=1,nx+1
            PTCec_old(i)=0.5d0*(PTCcc(i-1)+PTCcc(i))
            do n=1,Nspec
               rhoTDec_old(n,i)=0.5d0*(rhoTDcc(n,i-1)+rhoTDcc(n,i))
               rhoDiec_old(n,i)=0.5d0*(rhoDicc(n,i-1)+rhoDicc(n,i))
               do m=1,Nspec
                  rhoDijec_old(n,m,i)=0.5d0*(rhoDijcc(n,m,i-1)+rhoDijcc(n,m,i))
               enddo
            enddo
         enddo
         
         call LinOpApply(LofS,scal_old,PTCec_old,rhoTDec_old,rhoDijec_old,dx)

c     Form explicit update
         do i=1,nx+1
            do n=1,Nspec
               scal_new(FirstSpec+n-1,i)=scal_old(FirstSpec+n-1,i) +
     &              dt*LofS(n,i)
            enddo
            scal_new(RhoH,i)=scal_old(RhoH,i) + dt*LofS(Nspec+1,i)
            scal_new(Density,i) = scal_old(Density,i)
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
            maxIters = MAX(maxIters,Niter)
         enddo

c         print *,'t=',time,' RhoH->T maxIters: ',maxIters
         time = time + dt
         print *,'t=',time
      enddo
      call print_soln(time,scal_new,'soln_end.dat',dx,problo)

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
      integer i,n,m

      dxInv2 = 1.d0/(dx*dx)
      do i=0,nx+1
         do n=1,Nspec
            Y(n,i) = S(FirstSpec+n-1,i) / S(Density,i) 
         enddo
         call CKYTX(Y(1,i),IWRK,RWRK,X(1,i))
      enddo
      do i=1,nx+1
         do n=1,Nspec
            de(n,i) = (X(n,i)-X(n,i-1)) / dx
         enddo
         de(Nspec+1,i) = (S(Temp,i)-S(Temp,i-1)) / dx


         Te = 0.5d0*(S(Temp,i)+S(Temp,i-1))
         do n=1,Nspec
            Ye(n) = 0.5d0*(Y(n,i)+Y(n,i-1))
         enddo
         call CKMMWY(Ye,IWRK,RWRK,WWe)
         call CKHMS(Te,IWRK,RWRK,He)
         
         q(i) = 0.d0
         do n = 1,Nspec
            F(n,i) = - Ye(n)*rhoTDec(n,i)*de(Nspec+1,i)/Te
            do m = 1,Nspec
               F(n,i) = F(n,i) - rhoDijec(n,m,i)*de(m,i)
            enddo
            q(i) = q(i) - (RU*Te/WWe)*rhoTDec(n,i)*de(n,i) + He(n)*F(n,i)
         enddo
         q(i) = q(i) - PTCec(i)*de(Nspec+1,i)
      enddo

      do i=1,nx
         do n = 1,Nspec
            LofS(n,i) = - dxInv2*(F(n,i+1) - F(n,i))
         enddo
         LofS(Nspec+1,i) = - dxInv2*(q(i+1) - q(i))
      enddo

      end


      subroutine calc_beta(rho,T,rhoY,PTC,rhoTD,rhoDij,rhoDi)
      include 'spec.h'
      real*8 rho, T, rhoY(maxspec), Y(maxspec), X(maxspec), WW, CPMS,PTC
      real*8 rhoTD(maxspec),rhoDij(maxspec,maxspec),rhoDijt(maxspec*maxspec)
      real*8 rhoDi(maxspec)
      integer n,m,cnt

      do n=1,Nspec
         Y(n) = rhoY(n) / rho
      enddo
      call CKYTX(Y,IWRK,RWRK,X)
      call CKMMWY(Y,IWRK,RWRK,WW)
      call CKCPBS(T,Y,IWRK,RWRK,CPMS)
      call EGSPAR(T,X,Y,CPMS,EGRWRK,EGIWRK)
      call EGSLTDR5(T,Y,WW,EGRWRK,EGIWRK,PTC,rhoTD,rhoDijt)
      cnt = 1
      do n=1,Nspec
         do m=1,Nspec
            rhoDij(n,m) = rhoDijt(cnt)
            cnt = cnt+1
         enddo
      enddo
c     Mixture-averaged transport coefficients
      CALL EGSV1(Pcgs,T,Y,WW,EGRWRK,rhoDi)
      do n=1,Nspec
         rhoDi(n) = rho * WW * rhoDi(n) / mwt(n)
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

