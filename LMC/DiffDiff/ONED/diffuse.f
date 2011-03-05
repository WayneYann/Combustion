      subroutine mcdd_RhoH(S_new,S_old,dx,dt,theta,time,step)
      implicit none
      include 'spec.h'
      real*8 S_new(maxscal,0:nx+1)
      real*8 S_old(maxscal,0:nx+1)
      real*8 dx, dt, theta, time
      integer step
      real*8 S_star(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1)
      real*8 rhoDijec(maxspec,maxspec,1:nx+1)
      real*8 rhoDiec(maxspec,1:nx+1)
      real*8 cpicc(1:maxspec,0:nx+1)

      real*8 LofS_star(maxspec+1,1:nx)
      real*8 aofs(maxscal,1:nx)
      real*8 RHS(maxspec+1,1:nx)
      real*8 mass(maxspec), fac, dtDx2Inv, dt_temp
      real*8 err, L2err(maxscal), maxerr, prev, sum
      real*8 update(maxscal,1:nx)
      integer Niters, i, n, RhoH_to_Temp, iCorrect, idx, N1d
      real*8 a(N1dMAX), b(N1dMAX), c(N1dMAX), r(N1dMAX), v(N1dMAX), gam(N1dMAX)
      real*8 iterTol, Peos(1:nx), lambda, cpb
      real*8 tmp(1:nx), URfac, dtInv
      integer maxerrComp, firstPass
      parameter (iterTol=1.d-10)
      character*(50) itername, diffusename

      itername = 'iter'
      diffusename = 'diffuse'

c      Niters = RhoH_to_Temp(S_old)
c      if (Niters.lt.0) then
c         print *,'RhoH->Temp failed before predictor'
c         stop
c      endif

      dtInv = 1.d0 / dt
      do i=1,nx
         do n=1,maxscal
            aofs(n,i) = (S_old(n,i) - S_new(n,i)) * dtInv
         enddo
      enddo

      if ( (theta.lt.0.d0) .or. (theta.gt.1.d0) ) then
         print *,'diffuse_FULL: bad theta',theta
         stop
      endif

      if (theta.lt.1.d0) then
         call ecCoef_and_dt(S_old,PTCec,rhoTDec,rhoDijec,rhoDiec,cpicc,-1.d0,dx)
         call neg_divF_H(RHS,S_old,PTCec,rhoTDec,rhoDijec,dx)
         do i=1,nx
            do n=1,Nspec
               RHS(n,i) = RHS(n,i)*(1.d0-theta)*dt + S_old(FirstSpec+n-1,i)
     &               - aofs(FirstSpec+n-1,i)*dt
            enddo
            RHS(Nspec+1,i) = RHS(Nspec+1,i)*(1.d0-theta)*dt + S_old(RhoH,i)
     &           - aofs(RhoH,i)*dt
         enddo
      endif

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
         call neg_divF_H(LofS_star,S_star,PTCec,rhoTDec,rhoDijec,dx)

c     
c     R = b - A.S_new, then S_new = S_new + R/lambda
c
c     However, here, lambda comes from the mixture-averaged system:
c
c     d(rho.Yi)/dt = -Div(F)
c     d(rho.H)/dt  = -Div(Q)
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

            cpb = 0.d0
            do n=1,Nspec
               lambda=1.d0
               if (i.ne.1) then
                  lambda = lambda + dtDx2Inv*rhoDiec(n,i  )/S_star(Density,i)
               endif
               if (i.ne.nx) then
                  lambda = lambda + dtDx2Inv*rhoDiec(n,i+1)/S_star(Density,i)
               endif

               S_new(FirstSpec+n-1,i) = S_star(FirstSpec+n-1,i)
     &              + URfac*( RHS(n,i) - S_star(FirstSpec+n-1,i) + theta*dt*LofS_star(n,i) )/lambda

               cpb = cpb + cpicc(n,i)*S_star(FirstSpec+n-1,i)
            enddo
            cpb = cpb / S_star(Density,i)

            lambda = 1.d0
            if (i.ne.1) then
               lambda = lambda + dtDx2Inv*PTCec(i  )/(S_star(Density,i)*cpb)
            endif
            if (i.ne.nx) then
               lambda = lambda + dtDx2Inv*PTCec(i+1)/(S_star(Density,i)*cpb)
            endif

            S_new(RhoH,i) = S_star(RhoH,i)
     &           + URfac*( RHS(Nspec+1,i) - S_star(RhoH,i) + theta*dt*LofS_star(Nspec+1,i) )/lambda
         enddo

c     Diffuse state, compute errors/etc
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

            update(RhoH,i) = S_new(RhoH,i)-S_star(RhoH,i)
            err = ABS(update(RhoH,i))/typVal(RhoH)
            L2err(RhoH) = L2err(RhoH) + err*err

            Niters = RhoH_to_Temp(S_new)
            if (Niters.lt.0) then
               print *,'RhoH->Temp failed before predictor'
               stop
            endif

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

c         call print_soln(icorrect,DBLE(icorrect),S_new,itername,dx)
c         do i=1,nx
c            tmp(i) = S_new(Density,i)
c         enddo
c         call print_update(icorrect,DBLE(icorrect),diffuse,Peos,tmp,diffusename,dx)
      enddo
 100  end

      subroutine mcdd_Temp(S_new,S_old,dx,dt,theta,time,step)
      implicit none
      include 'spec.h'
      real*8 S_new(maxscal,0:nx+1)
      real*8 S_old(maxscal,0:nx+1)
      real*8 dx, dt, theta, time
      integer step
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
      character*(50) itername, diffusename

      itername = 'iter'
      diffusename = 'diffuse'

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
         print *,'diffuse_FULL1: bad theta',theta
         stop
      endif

      if (theta.lt.1.d0) then
         call ecCoef_and_dt(S_old,PTCec,rhoTDec,rhoDijec,rhoDiec,cpicc,-1.d0,dx)
         call neg_divF_T(RHS,S_old,PTCec,rhoTDec,rhoDijec,cpicc,dx)
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
         call neg_divF_T(LofS_star,S_star,PTCec,rhoTDec,rhoDijec,cpicc,dx)

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
c               lambda = lambda + PTCec(i  )*rhoCpInv
            endif
            if (i.ne.nx) then
               lambda = lambda + (PTCec(i+1) + RT)*rhoCpInv
c               lambda = lambda + PTCec(i+1)*rhoCpInv
            endif
            S_new(Temp,i) = S_star(Temp,i)
     &           + URfac*( RHS(Nspec+1,i) - S_star(Temp,i) + theta*dt*LofS_star(Nspec+1,i) )/lambda
         enddo

c     Diffuse state, compute errors/etc
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

c         call print_soln(icorrect,DBLE(icorrect),S_new,itername,dx)
c         do i=1,nx
c            tmp(i) = S_new(Density,i)
c         enddo
c         call print_update(icorrect,DBLE(icorrect),diffuse,Peos,tmp,diffusename,dx)
      enddo
 100  end


