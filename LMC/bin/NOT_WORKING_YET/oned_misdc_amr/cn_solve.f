      subroutine cn_solve(scal_new,alpha,beta_cc,Rhs,dx,dt,n,
     $                    be_cn_theta,rho_flag,is_vel,lo,hi,bc)
      implicit none
      include 'spec.h'

      integer rho_flag
      real*8 scal_new(-2:nfine+1,nscal)
      real*8    alpha( 0:nfine-1)
      real*8  beta_cc(-1:nfine  ,nscal)
      real*8      Rhs( 0:nfine-1)
      real*8 dx
      real*8 dt
      real*8 be_cn_theta
      logical is_vel
      integer lo,hi,bc(2)
      
      integer i,n
      real*8 a(nfine),b(nfine),c(nfine)
      real*8 r(nfine),u(nfine),gam(nfine)
      real*8 fac
      real*8 beta(0:nfine)
      integer n_solve

c     rho_flag used here to deal with rho scaling differences:
c     viscosity,conduction (rho_flag=1): rho.du/dt=Div(D.Grad(u))
c     mass (rho_flag=2): d(rho.u)/dt=Div(rho.D.Grad(u))

      if (coef_avg_harm.eq.1) then
         do i=lo,hi+1
            beta(i) = 2.0d0 / (1.d0/beta_cc(i,n)+1.d0/beta_cc(i-1,n))
         enddo
      else
         do i=lo,hi+1
            beta(i) = 0.5d0*(beta_cc(i,n) + beta_cc(i-1,n))
         enddo
      endif

      fac = be_cn_theta * dt / (dx*dx)

      do i=lo,hi
         u(i+1) = 0.d0
         r(i+1) = Rhs(i)
         a(i+1) = -fac*beta(i  )
         c(i+1) = -fac*beta(i+1)
         b(i+1) = alpha(i) - (a(i+1)+c(i+1))
c     dirichlet inflow using ghost cell value
         if (i.eq.lo) then
            a(i+1) = 0.d0
            if (rho_flag .eq. 2) then
               r(i+1) = r(i+1) + 
     $              fac*beta(i)*scal_new(-1,n)/scal_new(-1,Density)
            else
               r(i+1) = r(i+1) + 
     $              fac*beta(i)*scal_new(-1,n)
            end if
c     neumann outflow uses phi(hi) = phi(hi-1)
         else if (i.eq.hi) then
            c(i+1) = 0.d0
            b(i+1) = alpha(i) - (a(i+1)+c(i+1))
         end if


      enddo

      n_solve = hi-lo+1

      call tridiag(a,b,c,r,u,gam,n_solve)
      
      do i=lo,hi
         scal_new(i,n) = u(i+1)
      enddo
      
      if (rho_flag.eq.2) then
         do i=lo,hi
            scal_new(i,n) = scal_new(i,n) * scal_new(i,Density)
         enddo
      endif

      if (is_vel) then
         call set_bc_v(scal_new,lo,hi,bc)
      else
         call set_bc_s(scal_new,lo,hi,bc)
      end if
      
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

