      subroutine cn_solve(scal_new,alpha,beta_cc,Rhs,dx,dt,n,
     $                    be_cn_theta,rho_flag)
      implicit none
      include 'spec.h'

      integer rho_flag
      real*8 scal_new(-1:nx  ,*)
      real*8    alpha( 0:nx-1)
      real*8  beta_cc(-1:nx  ,*)
      real*8      Rhs( 0:nx-1)
      real*8 dx
      real*8 dt
      real*8 be_cn_theta
      
      integer i,n
      real*8 a(nx),b(nx),c(nx)
      real*8 r(nx),u(nx),gam(nx)
      real*8 fac
      real*8 beta(0:nx)

c     rho_flag used here to deal with rho scaling differences:
c     viscosity,conduction (rho_flag=1): rho.du/dt=Div(D.Grad(u))
c     mass (rho_flag=2): d(rho.u)/dt=Div(rho.D.Grad(u))

      if (coef_avg_harm.eq.1) then
         do i = 0,nx
            beta(i) = 2.0d0 / (1.d0/beta_cc(i,n)+1.d0/beta_cc(i-1,n))
         enddo
      else
         do i = 0,nx
            beta(i) = 0.5d0*(beta_cc(i,n) + beta_cc(i-1,n))
         enddo
      endif

      fac = be_cn_theta * dt / (dx*dx)

cc     homogeneous neumann inflow
c      do i = 0,nx-1
c         u(i+1) = 0.d0
c         r(i+1) = Rhs(i)
c         a(i+1) = -fac*beta(i  )
c         c(i+1) = -fac*beta(i+1)
c         if (i.eq.0   ) a(i+1) = 0.d0
c         if (i.eq.nx-1) c(i+1) =  0.d0
c         b(i+1) = alpha(i) - (a(i+1)+c(i+1))
c      enddo

c     dirichlet inflow using ghost cell value
      do i = 0,nx-1
         u(i+1) = 0.d0
         r(i+1) = Rhs(i)
         a(i+1) = -fac*beta(i  )
         c(i+1) = -fac*beta(i+1)
         if (i.eq.nx-1) c(i+1) =  0.d0
         b(i+1) = alpha(i) - (a(i+1)+c(i+1))
         if (i.eq.0   ) then
            a(i+1) = 0.d0
            if (rho_flag .eq. 2) then
               r(i+1) = r(i+1) + 
     $              fac*beta(i)*scal_new(-1,n)/scal_new(-1,Density)
            else
               r(i+1) = r(i+1) + 
     $              fac*beta(i)*scal_new(-1,n)
            end if
         end if
      enddo

      call tridiag(a,b,c,r,u,gam,nx)
      
      do i = 0,nx-1
         scal_new(i,n) = u(i+1)
      enddo
      
      if (rho_flag.eq.2) then
         do i = 0,nx-1
            scal_new(i,n) = scal_new(i,n) * scal_new(i,Density)
         enddo
      endif
      
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

      return
      end

