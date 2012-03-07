      subroutine probinit (domnlo, domnhi)
      implicit none

      include 'spec.h'
      double precision domnlo, domnhi
      integer nPMF, IWRK
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

      do n = 1,Nspec
         Y(n) = Y_bc(n,LorR)
      end do
      T = T_bc(LorR)
      rho = rho_bc(LorR)
      h = h_bc(LorR)
      u = u_bc(LorR)

      end

c----------------------------------------------------------------------

      subroutine initdata(vel,scal,I_R,dx)
      implicit none
      include 'spec.h'
      double precision    dx
      double precision   vel(-2:nx+1)
      double precision  scal(-2:nx+1,*)
      double precision   I_R(0:nx-1,0:*)

      double precision  x, rho, Y(maxspec), T, h
      double precision xPMFlo, xPMFhi
      double precision valsPMF(maxspec+3), RWRK, time, sum
      integer i, n, nPMF, IWRK

      write(*,*)'*** initdata *****'
      write(*,*)'Pcgs = ', Pcgs

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

      do i = 0,nx-1
         do n=0,Nspec
            I_R(i,n) = 0.d0
         enddo            
      enddo

      call set_bc_s(scal,dx,0.d0)
      call set_bc_v(vel,dx,0.d0)

      end


c----------------------------------------------------------------------

      subroutine set_bc_s(scal,dx,time)
      implicit none
      include 'spec.h'
      double precision     dx, time
      double precision   scal(-2:nx+1,*)
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
      double precision    vel(-2:nx+1)
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
