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
      if (LorR.eq.on_lo) then
         u = V_in
      else
         u = u_bc(LorR)
      endif
      end

c----------------------------------------------------------------------

      subroutine initdata(vel,scal,I_R,dx,lo,hi,bc)
      implicit none
      include 'spec.h'
      double precision  vel(0:nlevs-1,-2:nfine+1)
      double precision scal(0:nlevs-1,-2:nfine+1,nscal)
      double precision  I_R(0:nlevs-1,-1:nfine  ,0:Nspec)
      double precision   dx(0:nlevs-1)
      integer            lo(0:nlevs-1)
      integer            hi(0:nlevs-1)
      integer            bc(0:nlevs-1,2)

      double precision  x, rho, Y(Nspec), T, h
      double precision xPMFlo, xPMFhi
      double precision valsPMF(Nspec+3), RWRK, sum, sigma, cent
      integer i, n, l, nPMF, IWRK

      do l=0,nlevs-1

         do i=lo(l),hi(l)
            x = (dble(i)+0.5d0)*dx(l)

            if (probtype.eq.1) then
               xPMFlo = x - flame_offset - 0.5*dx(l)
               xPMFhi = x - flame_offset + 0.5*dx(l)
               call pmf(xPMFlo,xPMFhi,valsPMF,nPMF)
            else if (probtype.eq.2) then
               sigma = 0.0451d0
               cent = 0.d0
               valsPMF(1) = 5.d2 + 1.d3*dexp(-(x-cent)**2/(2.d0*sigma**2))
               valsPMF(2) = 0.d0
               do n=1,Nspec
                  valsPMF(3+n) = 0.d0
                  if (n.eq.iCH3OCH3) valsPMF(3+n) = .0374d0
                  if (n.eq.iCO2)     valsPMF(3+n) = .0522d0
                  if (n.eq.iO2)      valsPMF(3+n) = .1128d0
                  if (n.eq.iN2)      valsPMF(3+n) = .7192d0
                  if (n.eq.iH2O)     valsPMF(3+n) = .0784d0
               enddo
            endif

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
      call set_bc_s(scal(0,:,:),lo(0),hi(0),bc(0,:))
      call set_bc_v(vel(0,:)   ,lo(0),hi(0),bc(0,:))

      vel_TYP = ABS(vel(0,0))
      do l=0,nlevs-1
         do i=lo(l),hi(l)
            vel_TYP = MAX(vel_TYP,ABS(vel(l,i)))
         enddo
      enddo

      end

c----------------------------------------------------------------------

      subroutine set_bc_s(scal,lo,hi,bc)
      implicit none
      include 'spec.h'
      double precision scal(-2:nfine+1,nscal)
      integer lo,hi,bc(2)

!     local variables
      integer n, is, HorL
      double precision u, rho, Y(Nspec), T, hmix

      if (probtype.eq.1) then
         if (bc(1) .eq. 1) then
c     lo:  Dirichlet values for rho, Y, T, h
            HorL = on_lo
            call bcfunction(HorL, u, rho, Y, T, hmix)

            scal(lo-2:lo-1,Density) = rho
            scal(lo-2:lo-1,Temp) = T
            do n=1,Nspec
               is = FirstSpec + n - 1 
               scal(lo-2:lo-1,is) = Y(n) * rho
            enddo
            scal(lo-2:lo-1,RhoH) = hmix * rho
         end if
      else if (probtype.eq.2) then
c     lo:  Neumann for all 
         scal(lo-2:lo-1,Density) = scal(lo,Density)
         scal(lo-2:lo-1,Temp) = scal(lo,Temp)
         do n=1,Nspec
            is = FirstSpec + n - 1 
            scal(lo-2:lo-1,is) = scal(lo,is)
         enddo
         scal(lo-2:lo-1,RhoH) = scal(lo,RhoH)
      else
         print *,'Unknown probtype',probtype
         stop
      endif

      if (bc(2) .eq. 2) then
c     hi:  Neumann for all 
         scal(hi+1:hi+2,Density) = scal(hi,Density)
         scal(hi+1:hi+2,Temp) = scal(hi,Temp)
         do n=1,Nspec
            is = FirstSpec + n - 1 
            scal(hi+1:hi+2,is) = scal(hi,is)
         enddo
         scal(hi+1:hi+2,RhoH) = scal(hi,RhoH)
      end if

      end



c----------------------------------------------------------------------

      subroutine set_bc_v(vel,lo,hi,bc)
      implicit none
      include 'spec.h'
      double precision vel(-2:nfine+1)
      integer lo,hi,bc(2)

!     local variables
      integer HorL
      double precision u, rho, Y(Nspec), T, hmix

      if (bc(1) .eq. 1) then
c     lo:  Dirichlet
         HorL = on_lo
         call bcfunction(HorL, u, rho, Y, T, hmix)
         vel(lo-2:lo-1) = u
      end if

      if (bc(2) .eq. 2) then
c     hi:  Neumann for all
         vel(hi+1:hi+2) = vel(hi)
      end if

      end

c----------------------------------------------------------------------
