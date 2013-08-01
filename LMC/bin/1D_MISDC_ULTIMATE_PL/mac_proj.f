      subroutine macproj(macvel,rho,divu,dx,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8 macvel(0 :nfine)
      real*8    rho(-2:nfine+1)
      real*8   divu(-1:nfine)
      real*8 dx
      integer lo,hi,bc(2)

c     local variables
      integer i,n_solve

c     for tridiag solve
      real*8 a(nfine),b(nfine),c(nfine)
      real*8 r(nfine),u(nfine),gam(nfine)

      real*8 phi(-1:nfine+1)

      print *,'... mac_projection'

      do i=lo,hi
         u(i+1) = 0.d0
         r(i+1) = (macvel(i+1)-macvel(i))/dx - divu(i)
         a(i+1) = (2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)))
         b(i+1) = -(2.d0/dx**2)*
     &        (1.d0/(rho(i-1)+rho(i)) + 1.d0/(rho(i)+rho(i+1)))
         c(i+1) = (2.d0/dx**2)*(1.d0/(rho(i)+rho(i+1)))

         if (i .eq. lo .and. bc(1) .eq. 1) then
            a(i+1) = 0.d0
            b(i+1) = -(2.d0/dx**2)*(1.d0/(rho(i)+rho(i+1)))
            c(i+1) = -b(i+1)
         end if

         if (i .eq. hi .and. bc(2) .eq. 2) then
            a(i+1) = (1.d0/(3.d0*dx**2))*(1.d0/rho(i)) 
     &           + (2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)))
            b(i+1) = -(3.d0/dx**2)*(1.d0/rho(i))
     &           - (2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)))
            c(i+1) = 0.d0
         end if
      end do

      n_solve = hi-lo+1

      call tridiag(a,b,c,r,u,gam,n_solve)

      phi(lo:hi) = u(lo+1:hi+1)
      if (bc(1) .eq. 1) then
         phi(lo-1) = phi(lo)
      end if
      if (bc(2) .eq. 2) then
         phi(hi+1) = -2.d0*phi(hi) + (1.d0/3.d0)*phi(hi-1)
      end if

      do i=lo,hi+1
         macvel(i) = macvel(i) 
     &        - (2.d0/(rho(i-1)+rho(i)))*(phi(i)-phi(i-1))/dx
      end do
      
      end
