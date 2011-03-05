      

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

      integer Niter, maxIters, step, plot_int, do_advect
      real*8 res(NiterMAX)

c     Initialize chem/tran database
      call initchem()

c     Set defaults
      nsteps = 100
      plot_int = 100
      problo = 0.0d0
      probhi = 3.0d0
      flame_offset = 0.d0
      flame_offset = 1.d0

c      probhi = 1.d0
c      flame_offset = -.1d0

c     For LiDryer problems
c      problo = 2.0d0
c      probhi = 3.5d0
c      flame_offset = 0.d0

      Patm = 1.d0
      probtype = 1
      dtRedFac = 5.d0
      alt_spec_diffuse = 0
      mcdd_do_RhoH = 1
      do_advect = 1
      setTfromH = 1
      rhoInTrans = 1
      Ncorrect = 200
      Ncorrect = 1000
      be_cn_theta = 0.5d0
      outname = 'soln'

      call CKRP(IWRK,RWRK,RU,RUC,P1ATM)
      Pcgs = Patm * P1ATM

      Vel = 1
      Density = Vel + 1
      FirstSpec = Density + 1
      LastSpec = FirstSpec + Nspec - 1
      RhoH = LastSpec + 1
      Temp = RhoH + 1
      dx = (probhi-problo)/DBLE(nx)
      call init_soln(scal_new,time,dx)

      call print_soln(0,time,scal_new,outname,dx)

      dt = 1.d-5
      do step=1,nsteps
         do i=0,nx+1
            do n=1,maxscal
               scal_old(n,i) = scal_new(n,i)
            enddo
         enddo
         call apply_bcs(scal_old,time,step)

         if (do_advect.eq.1) then
            call advect(scal_new,scal_old,dx,dt,be_cn_theta,time,step)
         else
            do i=0,nx+1
               do n=1,maxscal
                  scal_new(n,i) = scal_old(n,i)
               enddo
            enddo
         endif

         if (mcdd_do_RhoH.eq.1) then
           call mcdd_RhoH(scal_new,scal_old,dx,dt,be_cn_theta,time,step)
         else
           call mcdd_Temp(scal_new,scal_old,dx,dt,be_cn_theta,time,step)
         endif

         time = time + dt
         print *,'step=', step, ' t=',time,' dt=',dt

         if (MOD(step,plot_int).eq.0  .or.  step.eq.nsteps) then
            call print_soln(step,time,scal_new,outname,dx)
         endif
      enddo
 100  continue
      end

      subroutine advect(S_new,S_old,dx,dt,theta,time,step)
      implicit none
      include 'spec.h'
      real*8 S_new(maxscal,0:nx+1)
      real*8 S_old(maxscal,0:nx+1)
      real*8 dx, dt, theta, time
      integer step, i, n, idx, Niters, RhoH_to_Temp
      real*8 SL(maxscal), SR(maxscal), aofs, dxInv

      if ( (theta.lt.0.d0) .or. (theta.gt.1.d0) ) then
         print *,'advect: bad theta',theta
         stop
      endif

      dxInv = 1.d0/dx
      do i=1,nx
c     Do simple upwind for now
         do n=1,maxscal
            SL(n) = S_old(n,i-1)
            SR(n) = S_old(n,i  )
         enddo

         S_new(Density,i) = 0.d0
         do n=1,Nspec
            idx = FirstSpec+n-1
            aofs = (SR(Vel)*SR(idx) - SL(Vel)*SL(idx)) * dxInv
            S_new(idx,i) = S_old(idx,i) - dt * aofs
            S_new(Density,i) = S_new(Density,i) + S_new(idx,i)
         enddo
         aofs = (SR(Vel)*SR(RhoH) - SL(Vel)*SL(RhoH)) * dxInv
         S_new(RhoH,i) = S_old(RhoH,i) - dt * aofs         
      enddo
      Niters = RhoH_to_Temp(S_old)
      if (Niters.lt.0) then
         print *,'RhoH->Temp failed in advect'
         stop
      endif
      end
