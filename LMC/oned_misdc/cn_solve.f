      subroutine cn_solve(nx,scal_new,alpha,beta_cc,Rhs,dx,dt,n,
     $                    be_cn_theta,rho_flag,rhohalf)
      implicit none
      include 'spec.h'

      integer nx, rho_flag
      real*8 scal_new(-1:nx  ,*)
      real*8    alpha( 0:nx-1)
      real*8  beta_cc(-1:nx  ,*)
      real*8      Rhs( 0:nx-1)
      real*8  rhohalf(0:*)
      real*8 dx
      real*8 dt
      real*8 be_cn_theta
      
      integer i,n
      real*8 a(nx),b(nx),c(nx)
      real*8 r(nx),u(nx)
      real*8 fac
      real*8 beta(0:nx)

c     rho_flag used here to deal with rho scaling differences:
c     viscosity,conduction (rho_flag=1): rho.du/dt=Div(D.Grad(u))
c     mass (rho_flag=2): d(rho.u)/dt=Div(rho.D.Grad(u))

      do i = 0,nx
         beta(i) = 2.0d0 /(1.d0/beta_cc(i,n)+1.d0/beta_cc(i-1,n))
      enddo

      fac = be_cn_theta * dt / (dx*dx)

      do i = 0,nx-1
         
         if (rho_flag.eq.3) then
            fac = fac / rhohalf(i)
         endif

         u(i+1) = 0.d0
         r(i+1) = Rhs(i)
         a(i+1) = -fac*beta(i  )
         c(i+1) = -fac*beta(i+1)
         if (i.eq.0   ) a(i+1) = 0.d0
         if (i.eq.nx-1) c(i+1) =  0.d0
         b(i+1) = alpha(i) - (a(i+1)+c(i+1))
      enddo

      call tridiag(a,b,c,r,u,nx)
      
      do i = 0,nx-1
         scal_new(i,n) = u(i+1)
      enddo
      
      do i = 0,nx-1
         if (rho_flag.eq.2) then
            scal_new(i,n) = scal_new(i,n) * scal_new(i,Density)
         endif
      enddo
      
      end

c *************************************************************************
c ** TRIDIAG **
c ** Do a tridiagonal solve 
c *************************************************************************

      subroutine tridiag(a,b,c,r,u,n)

      implicit none

      integer n
      real*8 a(n),b(n),c(n),r(n),u(n)
      integer j
      real*8 bet 
      real*8 gam(n)
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

