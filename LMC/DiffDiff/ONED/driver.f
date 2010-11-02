      

      subroutine lmc()
      implicit none
      include 'spec.h'
      integer nsteps

      real*8  scal_old(maxscal,0:nx+1)
      real*8  scal_new(maxscal,0:nx+1)
      real*8  LofS(maxspec+1,1:nx)

      real*8 dx, dt, enth, cpb
      real*8 rhom, rhop, rhoc, ym, yp, yc
      real*8 sum, maxsum, maxDiff, avgMag, T, rho, Yhalf
      integer i, Npmf, n, m
      real*8 x, time
      real*8 Patm, pmfdata(maxspec+3), mole(maxspec), mass(maxspec)
      real*8 rhoFactor
      character*(72) outname

      integer Niter, maxIters, step, plot_int
      real*8 res(NiterMAX)

c     Initialize chem/tran database
      call initchem()

c     Set defaults
      nsteps = 50
      nsteps = 2000
      plot_int = 10
      problo = 0.0d0
      probhi = 3.5d0
      flame_offset = 1.1d0
      Patm = 1.d0
      probtype = 2
      dtRedFac = 5.d0
      alt_spec_update = 0
      advance_RhoH = 1
      setTfromH = 2
      rhoInTrans = 1
      Ncorrect = 1000
      outname = 'soln'

      call CKRP(IWRK,RWRK,RU,RUC,P1ATM)
      Pcgs = Patm * P1ATM

      Density = 1
      FirstSpec = Density + 1
      LastSpec = FirstSpec + Nspec - 1
      RhoH = LastSpec + 1
      Temp = RhoH + 1
      dx = (probhi-problo)/nx
      call init_soln(scal_new,time,dx)

      call print_soln(0,time,scal_new,outname,dx,problo)

      dt = 1.d-4
      do step=1,nsteps
         do i=0,nx+1
            do n=1,maxscal
               scal_old(n,i) = scal_new(n,i)
            enddo
         enddo
         call apply_bcs(scal_old,time,step)

         if (advance_RhoH.eq.1) then
            call update_beImplicit(scal_new,scal_old,dx,dt,time,step)
c            call update_RhoH(scal_new,scal_old,dx,dt)
         else
            call update_beImplicit(scal_new,scal_old,dx,dt,time,step)
c            call update_Temp(scal_new,scal_old,dx,dt)
         endif

         time = time + dt
         print *,'step=', step, ' t=',time,' dt=',dt

         if (MOD(step,plot_int).eq.0  .or.  step.eq.nsteps) then
            call print_soln(step,time,scal_new,outname,dx,problo)
         endif
      enddo
 100  continue

      end

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

      subroutine LinOpApply(LofS,S,PTCec,rhoTDec,rhoDijec,dx)
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1),rhoDijec(maxspec,maxspec,1:nx+1),dx
      real*8 coef(maxspec+1,1:nx)
      real*8 Y(maxspec,0:nx+1), X(maxspec,0:nx+1), WWe, CPMS,PTC
      real*8 de(maxspec+1,1:nx+1), q(1:nx+1), F(maxspec,1:nx+1)
      real*8 Ye(maxspec), Te, dxInv2, He(maxspec)
      real*8 rhoe, enthe
      integer i,n,m,Niter,maxIter
      real*8 res(NiterMAX)

      dxInv2 = 1.d0/(dx*dx)
      do i=0,nx+1
         do n=1,Nspec
            Y(n,i) = S(FirstSpec+n-1,i) / S(Density,i) 
         enddo
         call CKYTX(Y(1,i),IWRK,RWRK,X(1,i))
      enddo

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

      subroutine LinOpApplyApprox(LofS,S,PTCec,rhoDiec,dx)
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoDiec(maxspec,1:nx+1),dx
      real*8 Y(maxspec,0:nx+1), X(maxspec,0:nx+1), WWe, CPBS,PTC
      real*8 de(maxspec+1,1:nx+1), q(1:nx+1), F(maxspec,1:nx+1)
      real*8 Ye(maxspec), Te, dxInv2, He(maxspec)
      real*8 rhoe, enthe
      integer i,n,m,Niter,maxIter
      real*8 res(NiterMAX), lamOverCp

      dxInv2 = 1.d0/(dx*dx)
      do i=0,nx+1
         do n=1,Nspec
            Y(n,i) = S(FirstSpec+n-1,i) / S(Density,i) 
         enddo
      enddo

      maxiter=0
      do i=1,nx+1
         do n=1,Nspec
            de(n,i) = (Y(n,i)-Y(n,i-1)) / dx
         enddo
         de(Nspec+1,i) = (S(RhoH,i)-S(RhoH,i-1)) / dx

         Te = 0.5d0*(S(Temp,i)+S(Temp,i-1))
         rhoe = 0.5d0*(S(Density,i)+S(Density,i-1))
         if (setTfromH.eq.1) then
            enthe = 0.5d0*(S(RhoH,i)+S(RhoH,i-1))/rhoe
            call FORT_TfromHYpt(Te,enthe,Ye,Nspec,errMax,NiterMAX,res,Niter)
            if (Niter.lt.0) then
               print *,'RhoH->T failed at i=',i,'in LinOpApplyApprox'
               stop
            endif
            maxiter=MAX(maxiter,Niter)
         endif
         call CKHMS(Te,IWRK,RWRK,He) 
         call CKCPBS(Te,Y(1,i),IWRK,RWRK,CPBS) 

         lamOverCp = PTCec(i) / CPBS

         q(i) = 0.d0
         do n = 1,Nspec
            F(n,i) = - rhoDiec(n,i)*de(n,i)
            q(i) = q(i) + He(n)*(rhoDiec(n,i) - lamOverCp)*de(n,i)
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

      subroutine LinOp1Apply(LofS,S,PTCec,rhoTDec,rhoDijec,cpb,dx)
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1),rhoDijec(maxspec,maxspec,1:nx+1),dx
      real*8 coef(maxspec+1,1:nx)
      real*8 Y(maxspec,0:nx+1), X(maxspec,0:nx+1), WWe, cpb(0:nx+1), PTC
      real*8 de(maxspec+1,1:nx+1), q(1:nx+1), F(maxspec,1:nx+1)
      real*8 Ye(maxspec), Te, dxInv2, cpi(maxspec)
      real*8 rhoe, enthe, Fnavg, gTavg
      integer i,n,m,Niter,maxIter
      real*8 res(NiterMAX)

      dxInv2 = 1.d0/(dx*dx)
      do i=0,nx+1
         do n=1,Nspec
            Y(n,i) = S(FirstSpec+n-1,i) / S(Density,i) 
         enddo
         call CKYTX(Y(1,i),IWRK,RWRK,X(1,i))
      enddo

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

         q(i) = 0.d0
         do n = 1,Nspec
            F(n,i) = - Ye(n)*rhoTDec(n,i)*de(Nspec+1,i)/Te
            do m = 1,Nspec
               F(n,i) = F(n,i) - rhoDijec(n,m,i)*de(m,i)
            enddo
            q(i) = q(i) - (RU*Te/WWe)*rhoTDec(n,i)*de(n,i)
         enddo
         q(i) = q(i) - PTCec(i)*de(Nspec+1,i)
      enddo

      do i=1,nx
         do n = 1,Nspec
            LofS(n,i) = - dxInv2*(F(n,i+1) - F(n,i))
         enddo
         call CKCPMS(S(Temp,i),IWRK,RWRK,cpi) 
         LofS(Nspec+1,i) = 0.d0
         do n = 1,Nspec
            Fnavg = 0.5d0*(F(n,i+1)+F(n,i))
            gTavg = 0.5d0*(de(Nspec+1,i+1)+de(Nspec+1,i))
            LofS(Nspec+1,i) = LofS(Nspec+1,i) - Fnavg*gTavg*cpi(n)
         enddo
         LofS(Nspec+1,i) = (LofS(Nspec+1,i) - dxInv2*(q(i+1) - q(i)))/(S(Density,i)*cpb(i))
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
         call EGSPAR(T,X,Y,CPMS,EGRWRK,EGIWRK)
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
      str_len = 1
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


      subroutine print_soln(step,time,scal,filename,dx,plo)
      include 'spec.h'
      integer step
      real*8 time, scal(maxscal,0:nx+1), dx, plo
      real*8 Peos(1:nx), Y(maxspec)
      character*(*) filename
      character*(100) fstr
      character*(100) fname
      character*(maxspnml) names(maxspec)
      integer n,i,j, get_spec_name, nlen(maxspec)
      do n=1,Nspec
         nlen = get_spec_name(names(n), n)
      enddo
      do i=1,nx
         do n=1,Nspec
            Y(n) = scal(FirstSpec+n-1,i)/scal(Density,i)
         enddo
         call CKPY(scal(Density,i),scal(Temp,i),Y,IWRK,RWRK,Peos(i))
      enddo
      write(fname,'(a,a,I0.6,a)') trim(filename),'_',step,'.dat'
      open(unit=12,file=trim(fname))
      write(fstr,'(a1,i2,a2)') '(',Nspec+2,'a)'
      write(12,fstr) 'VARIABLES=X Rho ',(trim(names(n)),' ',n=1,Nspec),
     &     ' RhoH Temp Peos'
      write(12,'(a,I0.6,a,I0.6,5(a,g20.12))') 'ZONE I=',nx,' T= "STEP=',step,
     &     ' time=',time,'" DATAPACKING=POINT STRANDID=1 SOLUTIONTIME=',
     &     time
      do i=1,nx
         write(12,'(50g20.12)') (i+0.5d0)*dx + plo, scal(1,i),
     &        (scal(1+n,i)/scal(1,i),n=1,Nspec),scal(Nspec+2,i),
     &        scal(Nspec+3,i),(Peos(i)-Pcgs)/Pcgs
      enddo
      close(12)
      end

      subroutine print_update(step,time,scal,Peos,rho,filename,dx,plo)
      include 'spec.h'
      integer step
      real*8 time, scal(maxscal,0:nx+1), Peos(1:nx), rho(1:nx), dx, plo
      character*(*) filename
      character*(100) fstr
      character*(100) fname
      character*(maxspnml) names(maxspec)
      integer n,i,j, get_spec_name, nlen(maxspec)
      do n=1,Nspec
         nlen = get_spec_name(names(n), n)
      enddo
      write(fname,'(a,a,I0.6,a)') trim(filename),'_',step,'.dat'
      open(unit=12,file=trim(fname))
      write(fstr,'(a1,i2,a2)') '(',Nspec+2,'a)'
      write(12,fstr) 'VARIABLES=X Rho ',(trim(names(n)),' ',n=1,Nspec),
     &     ' RhoH Temp Peos'
      write(12,'(a,I0.6,a,I0.6,5(a,g20.12))') 'ZONE I=',nx,' T= "STEP=',step,
     &     ' time=',time,'" DATAPACKING=POINT STRANDID=1 SOLUTIONTIME=',
     &     time
      do i=1,nx
         write(12,'(50g20.12)') (i+0.5d0)*dx + plo, scal(Density,i),
     &        (scal(FirstSpec+n-1,i)/rho(i),n=1,Nspec),scal(RhoH,i),
     &        scal(Nspec+3,i),(Peos(i)-Pcgs)/Pcgs
      enddo
      close(12)
      end

      subroutine init_soln(S,time,dx)
      implicit none
      include 'spec.h'

      real*8  S(maxscal,0:nx+1), time, dx, x
      integer i, n, Npmf
      real*8 pmfdata(maxspec+3), mole(maxspec), mass(maxspec)

      typVal(Temp)    = 0.d0
      typVal(RhoH)    = 0.d0
      typVal(Density) = 0.d0
      do i=1,nx
         if (probtype.eq.1) then
            x = (i+0.5d0)*dx - flame_offset
            call pmf(x,x,pmfdata,Npmf)
            if (Npmf.ne.Nspec+3) then
               print *,'mismatched pmf'
               stop
            endif
            S(Temp,i) = pmfdata(1)
            do n=1,Nspec
               mole(n) = pmfdata(3+n)
            enddo
         else
            x = (i+0.5d0)*dx
            S(Temp,i) = 298.d0
            do n=1,Nspec
               mole(n) = 0.d0
            enddo
            if (x.lt.0.5d0*(problo+probhi)) then
               mole(1) = 0.3d0
               mole(4) = 0.3d0
               mole(9) = 0.4d0
               S(Temp,i) = 298.d0
            else
               mole(4) = 0.21d0
               mole(9) = 0.79d0
               S(Temp,i) = 600.d0
               S(Temp,i) = 298.d0
            endif
         endif

         call CKXTY(mole,IWRK,RWRK,mass)         
         call CKRHOY(Pcgs,S(Temp,i),mass,IWRK,RWRK,S(Density,i))
         call CKHBMS(S(Temp,i),mass,IWRK,RWRK,S(RhoH,i))

         do n=1,Nspec
            S(FirstSpec+n-1,i) = mass(n) * S(Density,i)
         enddo
         S(RhoH,i) = S(RhoH,i) * S(Density,i)
         
         typVal(Temp)    = MAX(ABS(S(Temp,i)),   typVal(Temp))
         typVal(Density) = MAX(ABS(S(Density,i)),typVal(Density))
         typVal(RhoH   ) = MAX(ABS(S(RhoH,i)),   typVal(RhoH))
      enddo
      do n=1,Nspec
         typVal(FirstSpec+n-1) = typVal(Density)
      enddo
      time = 0.d0

      end

      subroutine apply_bcs(S,time,step)
      implicit none
      include 'spec.h'
      real*8 S(maxscal,0:nx+1)
      real*8 time
      integer step,n
c     Left boundary grow cell
      S(Density,0) = S(Density,1)
      S(Temp,0) = S(Temp,1)
      do n=1,Nspec
         S(FirstSpec+n-1,0) = S(FirstSpec+n-1,1)
      enddo
      S(RhoH,0) = S(RhoH,1)
      
c     Right boundary grow cell
      S(Density,nx+1) = S(Density,nx)
      S(Temp,nx+1) = S(Temp,nx)
      do n=1,Nspec
         S(FirstSpec+n-1,nx+1) = S(FirstSpec+n-1,nx)
      enddo
      S(RhoH,nx+1) = S(RhoH,nx)
      end

      subroutine update_RhoH(S_new,S_old,dx,dt)
      implicit none
      include 'spec.h'

      real*8  S_new(maxscal,0:nx+1)
      real*8  S_old(maxscal,0:nx+1)
      real*8  LofS(maxspec+1,1:nx)

      real*8  PTCec_old(1:nx+1)
      real*8  rhoTDec_old(maxspec,1:nx+1)
      real*8  rhoDijec_old(maxspec,maxspec,1:nx+1)
      real*8  rhoDiec_old(maxspec,1:nx+1)

      real*8 dx, dt, enth, cpicc(1:maxspec,0:nx+1)
      real*8 rhom, rhop, rhoc, ym, yp, yc
      real*8 sum, maxsum, maxDiff, avgMag, T, rho, Yhalf
      integer i, Npmf, n, m
      real*8 x, time, mass(maxspec)

      integer Niters, RhoH_to_Temp

c     positive dt signals that ecCoef_and_dt is to return stable dt
      dt = big
      call ecCoef_and_dt(S_old,PTCec_old,rhoTDec_old,rhoDijec_old,rhoDiec_old,cpicc,dt,dx)
      call LinOpApply(LofS,S_old,PTCec_old,rhoTDec_old,rhoDijec_old,dx)
      
c     Form explicit update
      do i=1,nx
         
         if (alt_spec_update .eq. 0) then
            
            do n=1,Nspec
               S_new(FirstSpec+n-1,i)=S_old(FirstSpec+n-1,i) +
     &              dt*LofS(n,i)
            enddo
            S_new(Density,i) = S_old(Density,i)
            
         else if (alt_spec_update .eq. 1) then
            
            sum = 0.d0
            do n=1,Nspec-1
               S_new(FirstSpec+n-1,i)=S_old(FirstSpec+n-1,i) +
     &              dt*LofS(n,i)
               sum = sum + S_new(FirstSpec+n-1,i)
            enddo
            S_new(Density,i) = S_old(Density,i)
            S_new(LastSpec,i) = S_new(Density,i) - sum
            
         else if (alt_spec_update .eq. 2) then
            
            sum = 0.d0
            do n=1,Nspec
               S_new(FirstSpec+n-1,i)=S_old(FirstSpec+n-1,i) +
     &              dt*LofS(n,i)
               sum = sum + S_new(FirstSpec+n-1,i)
            enddo
            S_new(Density,i) = sum
            
         else
            print *,'invalid value for alt_spec_update: ',alt_spec_update
         endif
         
         S_new(RhoH,i)=S_old(RhoH,i) + dt*LofS(Nspec+1,i)
      enddo
      
c     Recompute temperature
      Niters = RhoH_to_Temp(S_new)
      if (Niters.lt.0) then
         print *,'RhoH->Temp failed after explicit RhoH update'
         stop
      endif
      end

      integer function RhoH_to_Temp(S)
      implicit none
      include 'spec.h'
      real*8 S(maxscal,0:nx+1)
      real*8 mass(maxspec), enth
      integer Niter, maxIters, i, n
      real*8 res(NiterMAX)
      maxIters=0
      do i=1,nx
         do n=1,Nspec
            mass(n) = S(FirstSpec+n-1,i)/S(Density,i)
         enddo
         enth = S(RhoH,i)/S(Density,i)
         call FORT_TfromHYpt(S(Temp,i),enth,mass,
     &        Nspec,errMax,NiterMAX,res,Niter)
         if (Niter.lt.0) then
            print *,'RhoH->T failed at i=',i
            goto 100
         endif
         maxIters = MAX(maxIters,Niter)
      enddo
 100  RhoH_to_Temp = Niter
      end

      subroutine update_beImplicit(S_new,S_old,dx,dt,time,step)
      implicit none
      include 'spec.h'
      real*8 S_new(maxscal,0:nx+1)
      real*8 S_old(maxscal,0:nx+1)
      real*8 dx, dt, time, step
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1)
      real*8 rhoDijec(maxspec,maxspec,1:nx+1)
      real*8 rhoDiec(maxspec,1:nx+1)
      real*8 cpicc(1:maxspec,0:nx+1)

      real*8 LofS_old(maxspec+1,1:nx)
      real*8 LofS_new(maxspec+1,1:nx)

      real*8 mass(maxspec), fac, dtDxInv, dtDx2Inv, dt_temp, newVal
      real*8 err, L2err(maxscal), maxerr, prev, sum
      real*8 update(maxscal,1:nx), rho_new(1:nx)
      integer Niters, i, n, RhoH_to_Temp, iCorrect, idx, N1d
      real*8 a(N1dMAX), b(N1dMAX), c(N1dMAX), r(N1dMAX), v(N1dMAX), gam(N1dMAX)
      real*8 LT, RT, cpb, rhoCpInv, iterTol, Peos(1:nx)
      integer maxerrComp, relax_not_solve
      parameter (iterTol=1.d-10)

      character*(50) itername, updatename
      real*8 plo

      plo = 0.d0
      itername = 'iter'
      updatename = 'update'

      Niters = RhoH_to_Temp(S_old)
      if (Niters.lt.0) then
         print *,'RhoH->Temp failed before predictor'
         stop
      endif

      call apply_bcs(S_old,time,step)
      call ecCoef_and_dt(S_old,PTCec,rhoTDec,rhoDijec,rhoDiec,cpicc,-1.d0,dx)
      call LinOpApplyApprox(LofS_old,S_old,PTCec,rhoDiec,dx)

c     Initialize Snew to Sold
      do i=0,nx+1
         do n=1,maxscal
            S_new(n,i) = S_old(n,i)
         enddo
      enddo

c     Iteratively solve the following system
c
c     (rho.Y)_t = Div( rho.D.Grad(Yi) )
c
c     rho.cp.(T)_t = Div(lambda.Grad(T) ) + cpi.rho.Di.Grad(Yi).Grad(T)
c
c     Use Backward-Euler, tridiagonal solver, and 
c     (1) All gradients, transport coeffs evaluated on edges
c     (2) All cp, cpi evaluated on centers
c     (3) f=dt/dx2
c     (4) Lag coefficients, rho.Di.Grad(Yi) in T eqn, and cps
c

      call print_soln(0,0.d0,S_new,itername,dx,plo)

      relax_not_solve = 0

      do iCorrect=1,Ncorrect
            

         call apply_bcs(S_new,time,step)

         call ecCoef_and_dt(S_new,PTCec,rhoTDec,rhoDijec,rhoDiec,cpicc,-1.d0,dx)

         call build_approx_Ax_b(S_new,S_old,PTCec,rhoDiec,cpicc,dx,dt,time,step,
     &        a,b,c,r,N1d)
            
         if (relax_not_solve.eq.1) then

            do i=1,nx
               do n=1,Nspec
                  idx = nx*(n-1) + i
                  v(idx) = r(idx)
                  if (i.ne.1) then
                     v(idx) = v(idx) - a(idx)*S_new(FirstSpec+n-1,i-1)
                  endif
                  if (i.ne.nx) then
                     v(idx) = v(idx) - c(idx)*S_new(FirstSpec+n-1,i+1)
                  endif
                  v(idx) = v(idx) / b(idx)
               enddo
               idx = nx*Nspec + i
               v(idx) = r(idx)
               if (i.ne.1) then
                  v(idx) = v(idx) - a(idx)*S_new(Temp,i-1)
               endif
               if (i.ne.nx) then
                  v(idx) = v(idx) - c(idx)*S_new(Temp,i+1)
               endif
               v(idx) = v(idx) / b(idx)
            enddo

         else

            call tridiag(a,b,c,r,v,gam,N1d)

         endif

         do n=1,Nspec+3
            L2err(n) = 0.d0
         enddo

         do i=1,nx
            rho_new(i) = 0.d0
            do n=1,Nspec
               idx = nx*(n-1) + i
               rho_new(i) = rho_new(i) + v(idx)
            enddo
            do n=1,Nspec
               idx = nx*(n-1) + i
               prev = S_new(FirstSpec+n-1,i)
               S_new(FirstSpec+n-1,i) = v(idx)
               update(FirstSpec+n-1,i) = S_new(FirstSpec+n-1,i)-prev

               err = ABS(update(FirstSpec+n-1,i))/(typVal(FirstSpec+n-1)/typVal(Density))
               L2err(FirstSpec+n-1) = L2err(FirstSpec+n-1) + err*err

               mass(n) = S_new(FirstSpec+n-1,i) / rho_new(i)
            enddo
            idx = nx*Nspec + i
            prev = S_new(Temp,i)
            S_new(Temp,i) = v(idx)
            update(Temp,i) = S_new(Temp,i)-prev
            err = ABS(update(Temp,i))/typVal(Temp)
            L2err(Temp) = L2err(Temp) + err*err

            prev = S_new(RhoH,i)
            call CKHBMS(S_new(Temp,i),mass,IWRK,RWRK,S_new(RhoH,i))
            S_new(RhoH,i) = S_new(RhoH,i)*rho_new(i)
            update(RhoH,i) = S_new(RhoH,i)-prev
            err = ABS(update(RhoH,i))/typVal(RhoH)
            L2err(RhoH) = L2err(RhoH) + err*err

            prev = S_new(Density,i)
            S_new(Density,i) = rho_new(i)
            update(Density,i) = S_new(Density,i)-prev
            err = ABS(update(Density,i))/typVal(Density)
            L2err(Density) = L2err(Density) + err*err

            call CKPY(S_new(Density,i),S_new(Temp,i),mass,IWRK,RWRK,Peos(i))

         enddo

         maxerr = 0.d0
         maxerrComp = -1
         do n=1,Nspec+3
            L2err(n) = SQRT(L2err(n))
            if (L2err(n).gt.maxerr) then
               maxerr = L2err(n)
               maxerrComp = n
            endif
         enddo

         if (maxerr.lt.iterTol) then
            print *,'..converged after',iCorrect-1,'iters'
            goto 100
         else            
            print *, 'iter, err, component:',iCorrect,maxerr,maxerrComp
         endif


         call print_soln(icorrect,DBLE(icorrect),S_new,itername,dx,plo)
c         call print_update(icorrect,DBLE(icorrect),update,Peos,rho_new,updatename,dx,plo)

      enddo
 100  end


      subroutine build_approx_Ax_b(S_new,S_old,PTCec,rhoDiec,cpicc,dx,dt,time,step,
     &                             a,b,c,r,N1d)
      implicit none
      include 'spec.h'
      real*8 S_new(maxscal,0:nx+1)
      real*8 S_old(maxscal,0:nx+1)
      real*8 dx, dt, time
      real*8 PTCec(1:nx+1)
      real*8 rhoDiec(maxspec,1:nx+1)
      real*8 cpicc(1:maxspec,0:nx+1)
      real*8 a(N1dMAX), b(N1dMAX), c(N1dMAX), r(N1dMAX)
      integer step, N1d

      real*8 be_cn_theta, theta, dtDx2Inv, dt_temp
      integer Niters, i, n, RhoH_to_Temp, idx
      real*8 LT, RT, cpb, rhoCpInv, rhop, rhom

      be_cn_theta = 1.d0

      Niters = RhoH_to_Temp(S_old)
      if (Niters.lt.0) then
         print *,'RhoH->Temp failed before building matrix'
         stop
      endif

c     Build tridiagonal matrix form of the following system:
c
c     (rho.Y)_t = Div( rho.D.Grad(Yi) )
c
c     and
c
c     rho.cp.(T)_t = Div(lambda.Grad(T) ) + cpi.rho.Di.Grad(Yi).Grad(T)
c
c     Use Backward-Euler, tridiagonal solver, and 
c     (1) All gradients, transport coeffs evaluated on edges
c     (2) All cp, cpi evaluated on centers
c     (3) f=dt/dx2
c     (4) Lag coefficients, rho.Di.Grad(Yi) in T eqn, and cps
c
c     Requires that S_old, S_new have grow cells.  Upon entry, S_new
c     contains a guess for the solution state

      theta = be_cn_theta
      call apply_bcs(S_new,time,step)

      dt_temp = dt
      dtDx2Inv = dt_temp / (dx * dx)
      N1d = (Nspec+1)*nx

      do i=1,nx

         rhop = 0.5d0 * (S_new(Density,i+1)+S_new(Density,i))
         rhom = 0.5d0 * (S_new(Density,i-1)+S_new(Density,i))

         LT = 0.d0
         RT = 0.d0
         do n=1,Nspec
            
            idx = nx*(n-1) + i
            
            a(idx) = -dtDx2Inv*rhoDiec(n,i  )/rhop
            c(idx) = -dtDx2Inv*rhoDiec(n,i+1)/rhom
            r(idx) = S_old(FirstSpec+n-1,i)
            
            LT = LT + 
     &           0.25d0*rhoDiec(n,i  )*(cpicc(n,i-1)+cpicc(n,i  ))
     &           * ( S_new(FirstSpec+n-1,i  )/S_new(Density,i  )
     &           -   S_new(FirstSpec+n-1,i-1)/S_new(Density,i-1) )
            
            RT = RT +
     &           0.25d0*rhoDiec(n,i+1)*(cpicc(n,i  )+cpicc(n,i+1))
     &           * ( S_new(FirstSpec+n-1,i+1)/S_new(Density,i+1)
     &           -   S_new(FirstSpec+n-1,i  )/S_new(Density,i  ) )

         enddo

         idx = nx*Nspec + i
         
         cpb = 0.d0
         do n=1,Nspec
            cpb = cpb + cpicc(n,i)*S_new(FirstSpec+n-1,i)/S_new(Density,i)
         enddo
         rhoCpInv = 2.d0 / (cpb*S_old(Density,i))
         
         a(idx) = -dtDx2Inv*( PTCec(i  ) - LT)*rhoCpInv
         c(idx) = -dtDx2Inv*( PTCec(i+1) + RT)*rhoCpInv
         r(idx) = S_old(Temp,i)
         
         if (i.eq.1) then
            do n=1,Nspec+1
               idx = nx*(n-1) + i
               a(idx) = 0.d0
               b(idx) = 1.d0 - c(idx)
            enddo
         else if (i.eq.nx) then
            do n=1,Nspec+1
               idx = nx*(n-1) + i
               c(idx) = 0.d0
               b(idx) = 1.d0 - a(idx)
            enddo
         else
            do n=1,Nspec+1
               idx = nx*(n-1) + i
               b(idx) = 1.d0 - a(idx) - c(idx)
            enddo
         endif
      enddo
      end

c *************************************************************************
c ** TRIDIAG **
c ** Do a tridiagonal solve 
c *************************************************************************

      subroutine tridiag(a,b,c,r,u,gam,n)

      implicit none

      integer n
      real*8 a(n),b(n),c(n),r(n),u(n),gam(n)
      integer j
      real*8 bet 
      if (b(1) .eq. 0) print *,'CANT HAVE B(1) = ZERO'

      bet = b(1)
      u(1) = r(1)/bet

      do j = 2,n
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        if (bet .eq. 0) then
          print *,'TRIDIAG FAILED '
          stop
        endif
        u(j) = (r(j)-a(j)*u(j-1))/bet
      enddo

      do j = n-1,1,-1
        u(j) = u(j) - gam(j+1)*u(j+1)
      enddo

      end

      subroutine update_RhoH_pjImplicit(S_new,S_old,dx,dt,time,step)
      implicit none
      include 'spec.h'
      integer step
      real*8 S_new(maxscal,0:nx+1)
      real*8 S_old(maxscal,0:nx+1)
      real*8 dx, dt, time
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1)
      real*8 rhoDijec(maxspec,maxspec,1:nx+1)
      real*8 rhoDiec(maxspec,1:nx+1)
      real*8 LofS_old(maxspec+1,1:nx)
      real*8 LofS_new(maxspec+1,1:nx)
      real*8 cpb(0:nx+1)

      real*8 be_cn_theta, theta, mass(maxspec), fac, dtDxInv2, dt_temp, newVal
      real*8 err, L2err(maxscal),prev, Tprev(1:nx), URFac
      integer Niters, i, n, RhoH_to_Temp, iCorrect

      character*(50) junkname
      real*8 plo

      plo = 0.d0
      junkname = 'iter'

      be_cn_theta = 0.5d0
      URFac = 0.5d0

      Niters = RhoH_to_Temp(S_old)
      if (Niters.lt.0) then
         print *,'RhoH->Temp failed before predictor'
         stop
      endif

c     Initialize Snew to Sold
      do i=0,nx+1
         do n=1,maxscal
            S_new(n,i) = S_old(n,i)
         enddo
      enddo

      call ecCoef_and_dt(S_old,PTCec,rhoTDec,rhoDijec,rhoDiec,cpb,-1.d0,dx)
      call LinOpApply(LofS_old,S_old,PTCec,rhoTDec,rhoDijec,dx)

c     Relaxation: Jacobi iteration based on BE (theta=1) and CN (theta=0.5)
c
c                    Snew-Sold = dt*( (1-theta).Lo + theta*(L* - a.S* + a.Snew))
c
c     (1-dt.a.theta).Snew = Sold + dt.( (1-theta).Lo + theta.(L* - a.S*))
c
c     (1+fac).Snew = Sold + dt.( (1-theta).Lo + theta.L* ) + fac.S*,
c             where fac = -dt.a.theta = theta.(bL+bR).dt/dx2
c
c     Snew = (Sold + dt.( (1-theta).Lo + theta.L* ) + fac.S*))/(1+fac)
c
      call print_soln(0,0.d0,S_new,junkname,dx,plo)
      do iCorrect=1,Ncorrect
            
         theta = 1.d0
         URFac = 5.d-2

         call ecCoef_and_dt(S_new,PTCec,rhoTDec,rhoDijec,rhoDiec,cpb,-1.d0,dx)
         call LinOpApply(LofS_new,S_new,PTCec,rhoTDec,rhoDijec,dx)

         dt_temp = URFac*dt
         dtDxInv2 = dt_temp / (dx*dx)
         do n=1,Nspec+3
            L2err(n) = 0.d0
         enddo
         do i=1,nx
            do n=1,Nspec
               fac = (rhoDiec(n,i) + rhoDiec(n,i+1))*theta*dtDxInv2
               prev = S_new(FirstSpec+n-1,i)
               S_new(FirstSpec+n-1,i) = (S_old(FirstSpec+n-1,i) 
     &              + dt_temp*( (1.d0-theta)*LofS_old(n,i) + theta*LofS_new(n,i) )
     &              + fac*S_old(FirstSpec+n-1,i) ) / (1.d0 + fac)
               err = ABS(S_new(FirstSpec+n-1,i)-prev)/(typVal(FirstSpec+n-1)/typVal(Density))
               L2err(FirstSpec+n-1) = L2err(FirstSpec+n-1) + err*err
            enddo
            fac = (PTCec(i)/(S_old(Density,i)*cpb(i)) + PTCec(i+1)/(S_old(Density,i+1)*cpb(i+1)))*theta*dtDxInv2
            prev = S_new(RhoH,i)
            S_new(RhoH,i) = (S_old(RhoH,i) 
     &           + dt_temp*( (1.d0-theta)*LofS_old(Nspec+1,i) + theta*LofS_new(Nspec+1,i) )
     &           + fac*S_old(RhoH,i) ) / (1.d0 + fac)
            err = ABS(S_new(RhoH,i) - prev)/typVal(RhoH)
            L2err(RhoH) = L2Err(RhoH) + err*err
            S_new(Density,i) = S_old(Density,i)
         enddo
         do i=1,nx
            Tprev(i) = S_new(Temp,i)
         enddo
         call apply_bcs(S_new,time,step)

         Niters = RhoH_to_Temp(S_new)
         if (Niters.lt.0) then
            print *,'RhoH->Temp failed after corrector',icorrect
            stop
         endif
         do i=1,nx
            err = ABS(S_new(Temp,i) - Tprev(i))/typVal(Temp)
            L2err(Temp) = err*err
         enddo
         do n=1,Nspec+3
            L2err(n) = SQRT(L2err(n))
         enddo
         write(6,12) 'Err:',iCorrect,(L2err(n),n=1,Nspec+3)
 12      format(a,i3,12e9.2)

         call print_soln(icorrect,DBLE(icorrect),S_new,junkname,dx,plo)
      enddo
      end

      subroutine update_Temp(S_new,S_old,dx,dt)
      implicit none
      include 'spec.h'
      real*8 S_old(maxscal,0:nx+1)
      real*8 S_new(maxscal,0:nx+1)
      real*8 PTCec_old(1:nx+1)
      real*8 rhoTDec_old(maxspec,1:nx+1)
      real*8 rhoDijec_old(maxspec,maxspec,1:nx+1)
      real*8 rhoDiec_old(maxspec,1:nx+1)
      real*8 dx, dt

      real*8  LofS(maxspec+1,1:nx), mass(maxspec), sum, cpb(0:nx+1)
      integer i, n

c     positive dt signals that ecCoef_and_dt is to return stable dt
      dt = big
      call ecCoef_and_dt(S_old,PTCec_old,rhoTDec_old,rhoDijec_old,rhoDiec_old,cpb,dt,dx)
      call LinOp1Apply(LofS,S_old,PTCec_old,rhoTDec_old,rhoDijec_old,cpb,dx)
      
c     Form explicit update
      do i=1,nx
         
         if (alt_spec_update .eq. 0) then
            
            do n=1,Nspec
               S_new(FirstSpec+n-1,i)=S_old(FirstSpec+n-1,i) +
     &              dt*LofS(n,i)
            enddo
            S_new(Density,i) = S_old(Density,i)
            
         else if (alt_spec_update .eq. 1) then
            
            sum = 0.d0
            do n=1,Nspec-1
               S_new(FirstSpec+n-1,i)=S_old(FirstSpec+n-1,i) +
     &              dt*LofS(n,i)
               sum = sum + S_new(FirstSpec+n-1,i)
            enddo
            S_new(Density,i) = S_old(Density,i)
            S_new(LastSpec,i) = S_new(Density,i) - sum
            
         else if (alt_spec_update .eq. 2) then
            
            sum = 0.d0
            do n=1,Nspec
               S_new(FirstSpec+n-1,i)=S_old(FirstSpec+n-1,i) +
     &              dt*LofS(n,i)
               sum = sum + S_new(FirstSpec+n-1,i)
            enddo
            S_new(Density,i) = sum
            
         else
            print *,'invalid value for alt_spec_update: ',alt_spec_update
         endif
         
         S_new(Temp,i)=S_old(Temp,i) + dt*LofS(Nspec+1,i)
         do n=1,Nspec
            mass(n) = S_new(FirstSpec+n-1,i)/S_new(Density,i)
         enddo
         call CKHBMS(S_new(Temp,i),mass,IWRK,RWRK,S_new(RhoH,i))
         S_new(RhoH,i) = S_new(RhoH,i) * S_new(Density,i)
      enddo
      end



