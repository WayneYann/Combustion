      

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
      real*8 x, time, be_cn_theta
      real*8 Patm, pmfdata(maxspec+3), mole(maxspec), mass(maxspec)
      real*8 rhoFactor
      character*(72) outname

      integer Niter, maxIters, step, plot_int
      real*8 res(NiterMAX)

c     Initialize chem/tran database
      call initchem()

c     Set defaults
      nsteps = 1
      plot_int = 10
      problo = 0.0d0
      probhi = 3.5d0
      flame_offset = 1.1d0

c     For LiDryer problems
c      problo = 2.0d0
c      probhi = 3.5d0
c      flame_offset = 0.d0

      Patm = 1.d0
      probtype = 1
      dtRedFac = 5.d0
      alt_spec_update = 0
      advance_RhoH = 0
      setTfromH = 2
      rhoInTrans = 1
      Ncorrect = 200
      be_cn_theta = 0.5d0
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

      call print_soln(0,time,scal_new,outname,dx)

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
            call update_FULL(scal_new,scal_old,dx,dt,be_cn_theta,time,step)
c            call update_beFULL(scal_new,scal_old,dx,dt,time,step)
c            call update_beImplicit(scal_new,scal_old,dx,dt,time,step)
c            call update_Temp(scal_new,scal_old,dx,dt)
         endif

         time = time + dt
         print *,'step=', step, ' t=',time,' dt=',dt

         if (MOD(step,plot_int).eq.0  .or.  step.eq.nsteps) then
            call print_soln(step,time,scal_new,outname,dx)
         endif
      enddo
 100  continue

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

      subroutine LinOpApplyApprox(LofS,S,PTCec,rhoDiec,cpicc,dx)
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoDiec(maxspec,1:nx+1),cpicc(maxspec,0:nx+1),dx
      real*8 Y(maxspec,0:nx+1), X(maxspec,0:nx+1), WWe, PTC
      real*8 de(maxspec+1,1:nx+1), q(1:nx+1), F(maxspec,1:nx+1)
      real*8 Ye(maxspec), Te, dxInv2, He(maxspec)
      real*8 rhoe, enthe, cpb
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
         cpb = 0.d0
         do n=1,Nspec
            cpb = cpb + 0.5d0*(cpicc(n,i-1)+cpicc(n,i))
         enddo

         lamOverCp = PTCec(i) / cpb

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

      subroutine LinOp1Apply(LofS,S,PTCec,rhoTDec,rhoDijec,cpicc,dx)
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1),rhoDijec(maxspec,maxspec,1:nx+1),dx
      real*8 cpicc(maxspec,0:nx+1)
      real*8 Y(maxspec,0:nx+1), X(maxspec,0:nx+1), WWe, PTC
      real*8 de(maxspec+1,1:nx+1), q(1:nx+1), F(maxspec,1:nx+1)
      real*8 Ye(maxspec), Te, dxInv, sum
      real*8 rhoe, enthe, cpb
      integer i,n,m,Niter,maxIter
      real*8 res(NiterMAX)

      dxInv = 1.d0/dx
      do i=0,nx+1
         do n=1,Nspec
            Y(n,i) = S(FirstSpec+n-1,i) / S(Density,i) 
         enddo
         call CKYTX(Y(1,i),IWRK,RWRK,X(1,i))
      enddo

      maxiter=0
      do i=1,nx+1
         do n=1,Nspec
            de(n,i) = dxInv * (X(n,i)-X(n,i-1))
         enddo
         de(Nspec+1,i) = dxInv * (S(Temp,i)-S(Temp,i-1))

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
            LofS(n,i) = - dxInv*(F(n,i+1) - F(n,i))
         enddo
         sum = 0.d0
         cpb = 0.d0
         do n = 1,Nspec
            cpb = cpb + Y(n,i)*cpicc(n,i)
            sum = sum - 0.25d0*dxInv*(
     &           + F(n,i+1)*(cpicc(n,i  )+cpicc(n,i+1))*(S(Temp,i+1)-S(Temp,i  ) )
     &           + F(n,i  )*(cpicc(n,i-1)+cpicc(n,i  ))*(S(Temp,i  )-S(Temp,i-1) ) )
         enddo
         LofS(Nspec+1,i) = -(dxInv * (q(i+1) - q(i))  -  sum)/(S(Density,i) * cpb)
      enddo
      end

      subroutine LinOp1ApplyApprox(LofS,S,PTCec,rhoDiec,cpicc,dx)
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1), rhoDiec(maxspec,1:nx+1),cpicc(maxspec,0:nx+1),dx
      real*8 q(1:nx+1), F(maxspec,1:nx+1),dxInv,sum,cpb
      integer i,n,m
      dxInv = 1.d0/dx
      do i=1,nx+1
         do n = 1,Nspec
            F(n,i) = -rhoDiec(n,i)*dxInv
     &           * ( S(FirstSpec+n-1,i)/S(Density,i) - S(FirstSpec+n-1,i-1)/S(Density,i-1) )
         enddo
         q(i) = -PTCec(i) * dxInv * ( S(Temp,i) - S(Temp,i-1) )
      enddo

      do i=1,nx
         do n = 1,Nspec
            LofS(n,i) = -dxInv*(F(n,i+1) - F(n,i))
         enddo
         sum = 0.d0         
         cpb = 0.d0
         do n = 1,Nspec
            cpb = cpb + cpicc(n,i)*S(FirstSpec+n-1,i)/S(Density,i  )
            sum = sum - 0.25d0*dxInv*(
     &           + F(n,i+1)*(cpicc(n,i  )+cpicc(n,i+1))*(S(Temp,i+1)-S(Temp,i  ) )
     &           + F(n,i  )*(cpicc(n,i-1)+cpicc(n,i  ))*(S(Temp,i  )-S(Temp,i-1) ) )
         enddo
         LofS(Nspec+1,i) = -(dxInv * (q(i+1) - q(i))  -  sum)/(S(Density,i) * cpb)
      enddo

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

      real*8 LofS_new(maxspec+1,1:nx)
      real*8 mass(maxspec), fac, dtDx2Inv, dt_temp
      real*8 err, L2err(maxscal), maxerr, prev, sum
      real*8 update(maxscal,1:nx), rho_new(1:nx)
      integer Niters, i, n, RhoH_to_Temp, iCorrect, idx, N1d
      real*8 a(N1dMAX), b(N1dMAX), c(N1dMAX), r(N1dMAX), v(N1dMAX), gam(N1dMAX)
      real*8 iterTol, Peos(1:nx)
      integer maxerrComp, relax_not_solve, new_alg
      parameter (iterTol=1.d-10)

      character*(50) itername, updatename
      itername = 'iter'
      updatename = 'update'

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

c      call print_soln(0,0.d0,S_new,itername,dx)

      relax_not_solve = 1

      do iCorrect=1,Ncorrect
            
         call apply_bcs(S_new,time,step)

         call ecCoef_and_dt(S_new,PTCec,rhoTDec,rhoDijec,rhoDiec,cpicc,-1.d0,dx)

         call build_approx_Ax_b(S_new,S_old,PTCec,rhoDiec,cpicc,dx,dt,time,step,
     &        a,b,c,r,N1d)
            
         if (relax_not_solve.eq.1) then

            new_alg = 1

            if (new_alg .eq. 1) then

c     Compute R(S_new) = b - A.S_new
               do i=1,nx
                  do n=1,Nspec
                     idx = nx*(n-1) + i
                     gam(idx) = b(idx)*S_new(FirstSpec+n-1,i)
                     if (i.ne.1) then
                        gam(idx) = gam(idx) + a(idx)*S_new(FirstSpec+n-1,i-1)
                     endif
                     if (i.ne.nx) then
                        gam(idx) = gam(idx) + c(idx)*S_new(FirstSpec+n-1,i+1)
                     endif
                     gam(idx) = r(idx) - gam(idx)
                  enddo
                  idx = nx*Nspec + i
                  gam(idx) = b(idx)*S_new(Temp,i)
                  if (i.ne.1) then
                     gam(idx) = gam(idx) + a(idx)*S_new(Temp,i-1)
                  endif
                  if (i.ne.nx) then
                     gam(idx) = gam(idx) + c(idx)*S_new(Temp,i+1)
                  endif
                  gam(idx) = r(idx) - gam(idx)
               enddo

c     Now compute S_new <- S_new + R/lambda
               do i=1,nx
                  do n=1,Nspec
                     idx = nx*(n-1) + i
                     v(idx) = S_new(FirstSpec+n-1,i) + gam(idx)/b(idx)
                  enddo
                  idx = nx*Nspec + i
                  v(idx) = S_new(Temp,i) + gam(idx)/b(idx)
               enddo

            else

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

            endif

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


c         call print_soln(icorrect,DBLE(icorrect),S_new,itername,dx)
c         call print_update(icorrect,DBLE(icorrect),update,Peos,rho_new,updatename,dx)

      enddo
 100  end


      subroutine update_beFULL(S_new,S_old,dx,dt,time,step)
      implicit none
      include 'spec.h'
      real*8 S_new(maxscal,0:nx+1)
      real*8 S_old(maxscal,0:nx+1)
      real*8 dx, dt, time, step
      real*8 S_star(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1)
      real*8 rhoDijec(maxspec,maxspec,1:nx+1)
      real*8 rhoDiec(maxspec,1:nx+1)
      real*8 cpicc(1:maxspec,0:nx+1)

      real*8 LofS_star(maxspec+1,1:nx)
      real*8 mass(maxspec), fac, dtDx2Inv, dt_temp
      real*8 err, L2err(maxscal), maxerr, prev, sum
      real*8 update(maxscal,1:nx)
      integer Niters, i, n, RhoH_to_Temp, iCorrect, idx, N1d
      real*8 a(N1dMAX), b(N1dMAX), c(N1dMAX), r(N1dMAX), v(N1dMAX), gam(N1dMAX)
      real*8 iterTol, Peos(1:nx), cpnm, cpnp, lambda, LT, RT, cpb, rhoCpInv
      real*8 tmp(1:nx), URfac
      integer maxerrComp, firstPass
      parameter (iterTol=1.d-10)
      character*(50) itername, updatename

      itername = 'iter'
      updatename = 'update'

      Niters = RhoH_to_Temp(S_old)
      if (Niters.lt.0) then
         print *,'RhoH->Temp failed before predictor'
         stop
      endif

c     Initialize S_star to S_old
      do i=0,nx+1
         do n=1,maxscal
            S_new(n,i) = S_old(n,i)
         enddo
      enddo


c      do i=0,nx+1
c         tmp(0) = (DBLE(i)+0.5d0)*dx+problo
c         S_new(Temp,i) = 300.d0 + 1200.d0*(tmp(0)/(nx*dx))
c      enddo
c      call print_cp_prime(S_new,'cp.dat',dx)
c      stop

      call print_soln(0,0.d0,S_new,itername,dx)

      URfac = .75d0
      dtDx2Inv = dt/(dx*dx)
      firstPass = 1
      do iCorrect=1,Ncorrect
            
         do i=0,nx+1
            do n=1,maxscal
               S_star(n,i) = S_new(n,i)
            enddo
         enddo

         call apply_bcs(S_star,time,step)
         if (firstPass.eq.1) then
            call ecCoef_and_dt(S_star,PTCec,rhoTDec,rhoDijec,rhoDiec,cpicc,-1.d0,dx)
c            firstPass = 0
         endif
         call LinOp1Apply(LofS_star,S_star,PTCec,rhoTDec,rhoDijec,cpicc,dx)
c         call LinOp1ApplyApprox(LofS_star,S_star,PTCec,rhoDiec,cpicc,dx)

c     
c     R = b - A.S_new, then S_new = S_new + R/lambda
c
c     However, here, lambda comes from the mixture-averaged system:
c
c     d(rho.Yi)/dt = -Div(F)
c     rho.Cp.dT/dt = -Div(q) + rho.Di.Cpi.Grad(Yi).Grad(T)
c
c     The BE discretization is:
c
c               phi - dt.L(phi) = phi_old  (Ax=b)
c     
c     The point-jacobi iteration is phi <- phi^s + R^s / lambda^s
c
c     or, explicitly:   phi = phi^s + (phi_old - phi^s + dt*L(phi^s) / lambda
c
c     where lambda is the diagonal entry in the matrix A defined above
c
         do i=1,nx
c            rhom = 0.5d0*(S_star(Density,i-1)+S_star(Density,i  ))
c            rhop = 0.5d0*(S_star(Density,i  )+S_star(Density,i+1))

            LT = 0.d0
            RT = 0.d0

            do n=1,Nspec
               lambda=1.d0
               if (i.ne.1) then
                  cpnm = 0.5d0 * ( cpicc(n,i-1) + cpicc(n,i  ) )
                  lambda = lambda + dtDx2Inv*rhoDiec(n,i  )/S_star(Density,i)
c                  LT = LT + 0.5d0 * rhoDiec(n,i  ) * cpnm
c     &                 * ( S_star(FirstSpec+n-1,i) - S_star(FirstSpec+n-1,i-1) )/rhom
                  LT = LT + 0.5d0 * rhoDiec(n,i  ) * cpnm
     &                 * ( S_star(FirstSpec+n-1,i  )/S_star(Density,i)
     &                 -   S_star(FirstSpec+n-1,i-1)/S_star(Density,i-1) )
               endif
               if (i.ne.nx) then
                  cpnp = 0.5d0 * ( cpicc(n,i  ) + cpicc(n,i+1) )
                  lambda = lambda + dtDx2Inv*rhoDiec(n,i+1)/S_star(Density,i)
c                  RT = RT + 0.5d0 * rhoDiec(n,i+1) * cpnp
c     &                 * ( S_star(FirstSpec+n-1,i+1) - S_star(FirstSpec+n-1,i) )/rhop
                  RT = RT + 0.5d0 * rhoDiec(n,i+1) * cpnp
     &                 * ( S_star(FirstSpec+n-1,i+1)/S_star(Density,i+1)
     &                 -   S_star(FirstSpec+n-1,i  )/S_star(Density,i) )

               endif
               S_new(FirstSpec+n-1,i) = S_star(FirstSpec+n-1,i)
     &              + URfac*( S_old(FirstSpec+n-1,i) - S_star(FirstSpec+n-1,i) + dt*LofS_star(n,i) )/lambda
            enddo

            rhoCpInv = 0.d0
            do n=1,Nspec
               rhoCpInv = rhoCpInv + cpicc(n,i)*S_star(FirstSpec+n-1,i)
            enddo
            rhoCpInv = dtDx2Inv / rhoCpInv
            
            lambda = 1.d0
            if (i.ne.1) then
               lambda = lambda + (PTCec(i  ) - LT)*rhoCpInv
            endif
            if (i.ne.nx) then
               lambda = lambda + (PTCec(i+1) + RT)*rhoCpInv
            endif
            S_new(Temp,i) = S_star(Temp,i)
     &           + URfac*( S_old(Temp,i) - S_star(Temp,i) + dt*LofS_star(Nspec+1,i) )/lambda
         enddo

c     Update state, compute errors/etc
         do n=1,Nspec+3
            L2err(n) = 0.d0
         enddo
         
         do i=1,nx

            S_new(Density,i) = 0.d0
            do n=1,Nspec
               S_new(Density,i) = S_new(Density,i) + S_new(FirstSpec+n-1,i)
            enddo
            update(Density,i) = S_new(Density,i)-S_star(Density,i)
            err = ABS(update(Density,i))/typVal(Density)
            L2err(Density) = L2err(Density) + err*err

            do n=1,Nspec
               update(FirstSpec+n-1,i) = S_new(FirstSpec+n-1,i)-S_star(FirstSpec+n-1,i)
               err = ABS(update(FirstSpec+n-1,i))/(typVal(FirstSpec+n-1)/typVal(Density))
               L2err(FirstSpec+n-1) = L2err(FirstSpec+n-1) + err*err
               mass(n) = S_new(FirstSpec+n-1,i) / S_new(Density,i)
            enddo

            update(Temp,i) = S_new(Temp,i)-S_star(Temp,i)
            err = ABS(update(Temp,i))/typVal(Temp)
            L2err(Temp) = L2err(Temp) + err*err

            call CKHBMS(S_new(Temp,i),mass,IWRK,RWRK,S_new(RhoH,i))
            S_new(RhoH,i) = S_new(RhoH,i)*S_new(Density,i)
            update(RhoH,i) = S_new(RhoH,i)-S_star(RhoH,i)
            err = ABS(update(RhoH,i))/typVal(RhoH)
            L2err(RhoH) = L2err(RhoH) + err*err

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
            print *,'..converged after',iCorrect-1,
     &           'iters, (maxerr,iterTol): ',maxerr,iterTol
            goto 100
         else            
            print *, 'iter, err, component:',iCorrect,maxerr,maxerrComp
         endif

         call print_soln(icorrect,DBLE(icorrect),S_new,itername,dx)
         do i=1,nx
            tmp(i) = S_new(Density,i)
         enddo
         call print_update(icorrect,DBLE(icorrect),update,Peos,tmp,updatename,dx)
      enddo
 100  end


      subroutine update_FULL(S_new,S_old,dx,dt,theta,time,step)
      implicit none
      include 'spec.h'
      real*8 S_new(maxscal,0:nx+1)
      real*8 S_old(maxscal,0:nx+1)
      real*8 dx, dt, theta, time, step
      real*8 S_star(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1)
      real*8 rhoDijec(maxspec,maxspec,1:nx+1)
      real*8 rhoDiec(maxspec,1:nx+1)
      real*8 cpicc(1:maxspec,0:nx+1)

      real*8 LofS_star(maxspec+1,1:nx)
      real*8 RHS(maxspec+1,1:nx)
      real*8 mass(maxspec), fac, dtDx2Inv, dt_temp
      real*8 err, L2err(maxscal), maxerr, prev, sum
      real*8 update(maxscal,1:nx)
      integer Niters, i, n, RhoH_to_Temp, iCorrect, idx, N1d
      real*8 a(N1dMAX), b(N1dMAX), c(N1dMAX), r(N1dMAX), v(N1dMAX), gam(N1dMAX)
      real*8 iterTol, Peos(1:nx), cpnm, cpnp, lambda, LT, RT, cpb, rhoCpInv
      real*8 tmp(1:nx), URfac
      integer maxerrComp, firstPass
      parameter (iterTol=1.d-10)
      character*(50) itername, updatename

      itername = 'iter'
      updatename = 'update'

      Niters = RhoH_to_Temp(S_old)
      if (Niters.lt.0) then
         print *,'RhoH->Temp failed before predictor'
         stop
      endif

c     Initialize S_star to S_old
      do i=0,nx+1
         do n=1,maxscal
            S_new(n,i) = S_old(n,i)
         enddo
      enddo

      if ( (theta.lt.0.d0) .or. (theta.gt.1.d0) ) then
         print *,'update_FULL: bad theta',theta
         stop
      endif

      if (theta.lt.1.d0) then
         call ecCoef_and_dt(S_old,PTCec,rhoTDec,rhoDijec,rhoDiec,cpicc,-1.d0,dx)
         call LinOp1Apply(RHS,S_old,PTCec,rhoTDec,rhoDijec,cpicc,dx)
         do i=1,nx
            do n=1,Nspec
               RHS(n,i) = RHS(n,i)*(1.d0-theta)*dt + S_old(FirstSpec+n-1,i)
            enddo
            RHS(Nspec+1,i) = RHS(Nspec+1,i)*(1.d0-theta)*dt + S_old(Temp,i)
         enddo
      endif


c      do i=0,nx+1
c         tmp(0) = (DBLE(i)+0.5d0)*dx+problo
c         S_new(Temp,i) = 300.d0 + 1200.d0*(tmp(0)/(nx*dx))
c      enddo
c      call print_cp_prime(S_new,'cp.dat',dx)
c      stop

      call print_soln(0,0.d0,S_new,itername,dx)

      URfac = .75d0
      dtDx2Inv = theta*dt/(dx*dx)
      firstPass = 1
      do iCorrect=1,Ncorrect
            
         do i=0,nx+1
            do n=1,maxscal
               S_star(n,i) = S_new(n,i)
            enddo
         enddo

         call apply_bcs(S_star,time,step)
         if (firstPass.eq.1) then
            call ecCoef_and_dt(S_star,PTCec,rhoTDec,rhoDijec,rhoDiec,cpicc,-1.d0,dx)
c            firstPass = 0
         endif
         call LinOp1Apply(LofS_star,S_star,PTCec,rhoTDec,rhoDijec,cpicc,dx)
c         call LinOp1ApplyApprox(LofS_star,S_star,PTCec,rhoDiec,cpicc,dx)

c     
c     R = b - A.S_new, then S_new = S_new + R/lambda
c
c     However, here, lambda comes from the mixture-averaged system:
c
c     d(rho.Yi)/dt = -Div(F)
c     rho.Cp.dT/dt = -Div(q) + rho.Di.Cpi.Grad(Yi).Grad(T)
c
c     The CN discretization is:
c
c               phi - theta.dt.L(phi) = phi_old + (1-theta).dt.L(phi_old) = RHS (Ax=b)
c     
c     The point-jacobi iteration is phi <- phi^s + R^s / lambda^s
c
c     or, explicitly:   phi = phi^s + (RHS - phi^s + dt*L(phi^s) / lambda
c
c     where lambda is the diagonal entry in the matrix A defined above
c
         do i=1,nx

            LT = 0.d0
            RT = 0.d0

            do n=1,Nspec
               lambda=1.d0
               if (i.ne.1) then
                  cpnm = 0.5d0 * ( cpicc(n,i-1) + cpicc(n,i  ) )
                  lambda = lambda + dtDx2Inv*rhoDiec(n,i  )/S_star(Density,i)
                  LT = LT + 0.5d0 * rhoDiec(n,i  ) * cpnm
     &                 * ( S_star(FirstSpec+n-1,i  )/S_star(Density,i)
     &                 -   S_star(FirstSpec+n-1,i-1)/S_star(Density,i-1) )
               endif
               if (i.ne.nx) then
                  cpnp = 0.5d0 * ( cpicc(n,i  ) + cpicc(n,i+1) )
                  lambda = lambda + dtDx2Inv*rhoDiec(n,i+1)/S_star(Density,i)
                  RT = RT + 0.5d0 * rhoDiec(n,i+1) * cpnp
     &                 * ( S_star(FirstSpec+n-1,i+1)/S_star(Density,i+1)
     &                 -   S_star(FirstSpec+n-1,i  )/S_star(Density,i) )

               endif
               S_new(FirstSpec+n-1,i) = S_star(FirstSpec+n-1,i)
     &              + URfac*( RHS(n,i) - S_star(FirstSpec+n-1,i) + theta*dt*LofS_star(n,i) )/lambda
            enddo

            rhoCpInv = 0.d0
            do n=1,Nspec
               rhoCpInv = rhoCpInv + cpicc(n,i)*S_star(FirstSpec+n-1,i)
            enddo
            rhoCpInv = dtDx2Inv / rhoCpInv
            
            lambda = 1.d0
            if (i.ne.1) then
               lambda = lambda + (PTCec(i  ) - LT)*rhoCpInv
            endif
            if (i.ne.nx) then
               lambda = lambda + (PTCec(i+1) + RT)*rhoCpInv
            endif
            S_new(Temp,i) = S_star(Temp,i)
     &           + URfac*( RHS(Nspec+1,i) - S_star(Temp,i) + theta*dt*LofS_star(Nspec+1,i) )/lambda
         enddo

c     Update state, compute errors/etc
         do n=1,Nspec+3
            L2err(n) = 0.d0
         enddo
         
         do i=1,nx

            S_new(Density,i) = 0.d0
            do n=1,Nspec
               S_new(Density,i) = S_new(Density,i) + S_new(FirstSpec+n-1,i)
            enddo
            update(Density,i) = S_new(Density,i)-S_star(Density,i)
            err = ABS(update(Density,i))/typVal(Density)
            L2err(Density) = L2err(Density) + err*err

            do n=1,Nspec
               update(FirstSpec+n-1,i) = S_new(FirstSpec+n-1,i)-S_star(FirstSpec+n-1,i)
               err = ABS(update(FirstSpec+n-1,i))/(typVal(FirstSpec+n-1)/typVal(Density))
               L2err(FirstSpec+n-1) = L2err(FirstSpec+n-1) + err*err
               mass(n) = S_new(FirstSpec+n-1,i) / S_new(Density,i)
            enddo

            update(Temp,i) = S_new(Temp,i)-S_star(Temp,i)
            err = ABS(update(Temp,i))/typVal(Temp)
            L2err(Temp) = L2err(Temp) + err*err

            call CKHBMS(S_new(Temp,i),mass,IWRK,RWRK,S_new(RhoH,i))
            S_new(RhoH,i) = S_new(RhoH,i)*S_new(Density,i)
            update(RhoH,i) = S_new(RhoH,i)-S_star(RhoH,i)
            err = ABS(update(RhoH,i))/typVal(RhoH)
            L2err(RhoH) = L2err(RhoH) + err*err

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
            print *,'..converged after',iCorrect-1,
     &           'iters, (maxerr,iterTol): ',maxerr,iterTol
            goto 100
         else            
            print *, 'iter, err, component:',iCorrect,maxerr,maxerrComp
         endif

         call print_soln(icorrect,DBLE(icorrect),S_new,itername,dx)
         do i=1,nx
            tmp(i) = S_new(Density,i)
         enddo
         call print_update(icorrect,DBLE(icorrect),update,Peos,tmp,updatename,dx)
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
            
            a(idx) = -dtDx2Inv*rhoDiec(n,i  )/rhom
            c(idx) = -dtDx2Inv*rhoDiec(n,i+1)/rhop
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
         rhoCpInv = 1.d0 / (cpb*S_new(Density,i))
         
         a(idx) = - dtDx2Inv*(PTCec(i  ) - LT)*rhoCpInv
         c(idx) = - dtDx2Inv*(PTCec(i+1) + RT)*rhoCpInv
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
      call print_soln(0,0.d0,S_new,junkname,dx)
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

         call print_soln(icorrect,DBLE(icorrect),S_new,junkname,dx)
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



