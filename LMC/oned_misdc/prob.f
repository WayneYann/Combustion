      subroutine probinit (domnlo, domnhi)
      implicit none

      include 'spec.h'
      double precision domnlo, domnhi
      integer untin, nPMF, IWRK
      double precision xPMF, valsPMF(maxspec+3), RWRK

c     Set ordering of variables in state
      Density = 1
      Temp = 2
      RhoH = 3
      RhoRT = 4
      FirstSpec = 5
      LastSpec = FirstSpec + (nspec-1)
      nscal = LastSpec

c     Set boundary data
      xPMF = domnlo - flame_offset
      call pmf(xPMF,xPMF,valsPMF,nPMF)
      call CKXTY(valsPMF(4),IWRK,RWRK,Y_bc(1,0))
      T_bc(0) = valsPMF(1)
      if (ABS(V_in) .lt. 1.d20) then
         u_bc(0) = V_in
      else
         u_bc(0) = valsPMF(2)
      endif

      call CKHBMS(T_bc(0),Y_bc(1,0),IWRK,RWRK,h_bc(0))
      call CKRHOY(Pcgs,T_bc(0),Y_bc(1,0),IWRK,RWRK,rho_bc(0))
      end

c-----------------------------------------------------------------------

      subroutine bcfunction(LorR, time, u, rho, Y, T, h)
      implicit none
      include 'spec.h'
      integer LorR, n
      double precision time, u, rho, Y(*), T, h

      if ( probtype.eq.1 ) then

         do n = 1,Nspec
            Y(n) = Y_bc(n,LorR)
         end do
         T = T_bc(LorR)
         rho = rho_bc(LorR)
         h = h_bc(LorR)
         u = u_bc(LorR)

      else
         print *,'bcfunction: invalid probtype'
         stop
      end if
      end

c----------------------------------------------------------------------

      subroutine initdata(vel,scal,I_R,dx)
      implicit none
      include 'spec.h'
      double precision    dx
      double precision   vel(-1:nx)
      double precision  scal(-1:nx  ,*)
      double precision   I_R(0:nx-1,0:*)

      double precision  x, rho, Y(maxspec), T, h
      double precision xPMFlo, xPMFhi
      double precision Z(0:maxspec),ZP(0:maxspec)
      double precision valsPMF(maxspec+3), RWRK, time, sum
      integer i, n, nPMF, IWRK

      write(*,*)'*** initdata *****'
      write(*,*)'Pcgs = ', Pcgs

      if (probtype.eq.1) then

         time = 0.d0
         do i = 0,nx-1
            x = (DBLE(i)+0.5D0)*dx
            
            xPMFlo = x - flame_offset - 0.5*dx
            xPMFhi = x - flame_offset + 0.5*dx
            call pmf(xPMFlo,xPMFhi,valsPMF,nPMF)

            call CKXTY(valsPMF(4),IWRK,RWRK,Y)
            T = valsPMF(1)

            if (iN2.gt.0  .and.  iN2.le.Nspec) then
               sum = 0.d0
               do n=1,Nspec
                  sum = sum + Y(n)
               enddo
               Y(iN2) = Y(iN2) - 1.d0 + sum
            endif

            call CKHBMS(T,Y,IWRK,RWRK,h)
            call CKRHOY(Pcgs,T,Y,IWRK,RWRK,rho)
            do n=1,Nspec
               scal(i,FirstSpec+n-1) = rho * Y(n)
            enddo
            scal(i,Density) = rho
            scal(i,Temp) = T
            scal(i,RhoH) = rho * h
            vel(i) = valsPMF(2)
            if (i.eq.0) then
               hmix_TYP = ABS(h)
            else
               hmix_TYP = MAX(hmix_TYP,ABS(h))
            endif
         enddo

         do n = 0,Nspec
            c_0(n) = 0.d0
            c_1(n) = 0.d0
         enddo
         if (nochem_hack .or. use_strang) then
            do i = 0,nx-1
               do n=0,Nspec
                  I_R(i,n) = 0.d0
               enddo            
            enddo
         else
            do i = 0,nx-1
               Z(0) = scal(i,Temp)
               do n=1,Nspec
                  Z(n) = scal(i,FirstSpec+n-1)
               enddo
               rhoh_INIT = scal(i,RhoH)
               call vodeF_T_RhoY(Nspec+1,time,Z(0),ZP(0),RWRK,IWRK)
               do n=0,Nspec
                  I_R(i,n) = ZP(N)
               enddo
            enddo
         endif
      else
         print *,'bcfunction: invalid probtype'
         stop
      endif
      end


c----------------------------------------------------------------------

      subroutine set_bc_s(scal,dx,time)
      implicit none
      include 'spec.h'
      double precision     dx, time
      double precision   scal(-1:nx  ,*)
      integer n, is, HorL
      double precision u, rho, Y(maxspec), T, hmix

c     Sets the grow cell to hold the boundary value
c     For Dirichlet condition, this value is to be applied at the cell

c     lo:  Dirichlet values for u, T, Y, compute rho, h
      HorL = on_lo
      call bcfunction(HorL, time, u, rho, Y, T, hmix)
      scal(-1,Density) = rho
      scal(-1,Temp) = T
      do n=1,Nspec
         is = FirstSpec + n - 1 
         scal(-1,is) = Y(n) * rho
      enddo
      scal(-1,RhoH) = hmix * rho
      
c     hi:  Neumann for all 
      scal(nx,Density) = scal(nx-1,Density)
      scal(nx,Temp) = scal(nx-1,Temp)
      do n=1,Nspec
         is = FirstSpec + n - 1 
         scal(nx,is) = scal(nx-1,is)
      enddo
      scal(nx,RhoH) = scal(nx-1,RhoH)
      end


c----------------------------------------------------------------------

      subroutine set_bc_v(vel,dx,time)
      implicit none
      include 'spec.h'
      integer              HorL
      double precision     dx, time
      double precision    vel(-1:nx)
      double precision u, rho, Y(maxspec), T, hmix

c     Sets the grow cell to hold the boundary value
c     For Dirichlet condition, this value is to be applied at the cell

c     lo:  Dirichlet values for u, T, Y, compute rho, h
c     hi:  Neumann

      HorL = on_lo
      call bcfunction(HorL, time, u, rho, Y, T, hmix)
      vel(-1) = u
      vel(nx) = vel(nx-1)
      end


c----------------------------------------------------------------------

      subroutine set_bc_grow_s(scal,dx,time)
      implicit none
      include 'spec.h'
      double precision     dx, time
      double precision   scal(-1:nx  ,*)
      integer n, i, is, HorL, leni, ib, IWRK
      double precision x(-1:max_order-2), xInt, RWRK
      double precision coef(-1:max_order-2)
      double precision Text, Yext(maxspec), hext, rext, rho_i

c     Set coeffs for polynomial extrap to fill grow cells
      leni = max_order-2
      do i=0,leni
         x(i) = i + 0.5D0
      end do
      x(-1) = 0.d0
      xInt = - 0.5d0
      call polyInterpCoeff(xInt, x, max_order, coef)

c     Set Dirichlet values into grow cells
      call set_bc_s(scal,dx,time)

c     Do extrap
      ib = -1
      do n=1,Nspec
         is = FirstSpec + n - 1
         Yext(n) = coef(-1) * scal(ib,is)/scal(ib,Density)
      enddo
      Text = coef(-1) * scal(ib,Temp)
      do i=0,leni
         rho_i = 0.d0
         do n=1,Nspec
            is = FirstSpec + n - 1
            rho_i = rho_i + scal(i,is)
         enddo
         do n=1,Nspec
            is = FirstSpec + n - 1
            Yext(n) = Yext(n) + scal(i,is)*coef(i)/rho_i
         enddo
         Text = Text + scal(i,Temp)*coef(i)
      enddo
      
c     On lo boundary grow cell, now have Y and T in state, get rho and h and fix
      call CKHBMS(Text,Yext,IWRK,RWRK,hext)
      call CKRHOY(Pcgs,Text,Yext,IWRK,RWRK,rext)
      
      scal(-1,Density) = rext
      do n=1,Nspec
         is = FirstSpec + n - 1
         scal(-1,is) = Yext(n) * rext
      enddo
      scal(-1,RhoH) = hext * rext
      scal(-1,Temp) = Text         

      end

c----------------------------------------------------------------------

      subroutine set_bc_grow_v(vel,dx,time)
      implicit none
      include 'spec.h'
      double precision     dx, time
      double precision    vel(-1:nx)
      integer n, i, is, HorL, leni, ib, IWRK
      double precision x(-1:max_order-2), xInt
      double precision coef(-1:max_order-2), vext


c     Set coeffs for polynomial extrap to fill grow cells
      leni = max_order-2
      do i=0,leni
         x(i) = i + 0.5D0
      end do
      x(-1) = 0.d0
      xInt = - 0.5d0
      call polyInterpCoeff(xInt, x, max_order, coef)

c     Set Dirichlet values into grow cells
      call set_bc_v(vel,dx,time)

c     Do extrap
      ib = -1
      vext = coef(-1) * vel(ib)
      do i=0,leni
         vext = vext + vel(i)*coef(i)
      enddo
      vel(-1) = vext
      end

c----------------------------------------------------------------------
c
c     polyInterpCoeff:
c  
c     This routine returns the Lagrange interpolating coefficients for a
c     polynomial through N points, evaluated at xInt (see Numerical Recipes,
c     v2, p102, e.g.):
c
c            (x-x2)(x-x3)...(x-xN)              (x-x1)(x-x2)...(x-x(N-1))
c    P(x) = ----------------------- y1  + ... + ------------------------  yN
c           (x1-x2)(x1-x3)...(x1-xN)            (x1-x2)(x1-x3)...(x1-xN)
c
c     P(xInt) = sum_(i=1)^(N) y[i]*c[i]
c
      subroutine polyInterpCoeff(xInt, x, N, c)
      implicit none
      integer N, i, j
      double precision xInt, x(N), c(N), num, den
      do j=1,N
         num = 1.d0
         den = 1.d0
         do i = 1,j-1
            num = num*(xInt - x(i))
            den = den*(x(j) - x(i))
         end do
         do i = j+1,N
            num = num*(xInt - x(i))
            den = den*(x(j) - x(i))
         end do
         c(j) = num/den
      end do
      end
