      subroutine probinit (domnlo)
      implicit none

      include 'spec.h'
      double precision domnlo
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

      subroutine bcfunction(LorR, u, rho, Y, T, h)
      implicit none
      include 'spec.h'
      integer LorR, n
      double precision u, rho, Y(*), T, h

      do n = 1,Nspec
         Y(n) = Y_bc(n,LorR)
      end do
      T = T_bc(LorR)
      rho = rho_bc(LorR)
      h = h_bc(LorR)
      u = u_bc(LorR)

      end

c----------------------------------------------------------------------

      subroutine initdata(vel,scal,I_R,dx,lo,hi)
      implicit none
      include 'spec.h'
      double precision  vel(0:nlevs-1,-2:nfine+1)
      double precision scal(0:nlevs-1,-2:nfine+1,nscal)
      double precision  I_R(0:nlevs-1,-1:nfine  ,0:Nspec)
      double precision   dx(0:nlevs-1)
      integer            lo(0:nlevs-1)
      integer            hi(0:nlevs-1)

      double precision  x, rho, Y(Nspec), T, h
      double precision xPMFlo, xPMFhi
      double precision valsPMF(Nspec+3), RWRK, sum
      integer i, n, l, nPMF, IWRK

      do l=0,nlevs-1

         do i=lo(l),hi(l)
            x = (dble(i)+0.5d0)*dx(l)

            xPMFlo = x - flame_offset - 0.5*dx(l)
            xPMFhi = x - flame_offset + 0.5*dx(l)
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
               scal(l,i,FirstSpec+n-1) = rho * Y(n)
            enddo
            scal(l,i,Density) = rho
            scal(l,i,Temp) = T
            scal(l,i,RhoH) = rho * h
            vel(l,i) = valsPMF(2)
            if (l .eq. 0 .and. i .eq. lo(0)) then
               hmix_TYP = ABS(h)
            else
               hmix_TYP = MAX(hmix_TYP,ABS(h))
            endif
         enddo

      enddo

      c_0 = 0.d0
      c_1 = 0.d0

      I_R = 0.d0

!     fill coarse level ghost cells
      call set_bc_s(scal(0,:,:))
      call set_bc_v(vel(0,:))

      vel_TYP = ABS(vel(0,0))
      do l=0,nlevs-1
         do i=lo(l),hi(l)
            vel_TYP = MAX(vel_TYP,ABS(vel(l,i)))
         enddo
      enddo

      end

c----------------------------------------------------------------------

      subroutine set_bc_s(scal)
      implicit none
      include 'spec.h'
      double precision scal(-2:nfine+1,nscal)

!     local variables
      integer n, is, HorL
      double precision u, rho, Y(Nspec), T, hmix

c     lo:  Dirichlet values for rho, Y, T, h
      HorL = on_lo
      call bcfunction(HorL, u, rho, Y, T, hmix)
      scal(-2:-1,Density) = rho
      scal(-2:-1,Temp) = T
      do n=1,Nspec
         is = FirstSpec + n - 1 
         scal(-2:-1,is) = Y(n) * rho
      enddo
      scal(-2:-1,RhoH) = hmix * rho
      
c     hi:  Neumann for all 
      scal(nx:nx+1,Density) = scal(nx-1,Density)
      scal(nx:nx+1,Temp) = scal(nx-1,Temp)
      do n=1,Nspec
         is = FirstSpec + n - 1 
         scal(nx:nx+1,is) = scal(nx-1,is)
      enddo
      scal(nx:nx+1,RhoH) = scal(nx-1,RhoH)
      end


c----------------------------------------------------------------------

      subroutine set_bc_v(vel)
      implicit none
      include 'spec.h'
      double precision vel(-2:nfine+1)

!     local variables
      integer HorL
      double precision u, rho, Y(Nspec), T, hmix

c     lo:  Dirichlet
      HorL = on_lo
      call bcfunction(HorL, u, rho, Y, T, hmix)
      vel(-2:-1) = u

c     hi:  Neumann for all
      vel(nx:nx+1) = vel(nx-1)
      end

c----------------------------------------------------------------------
