      subroutine update_rho(nx,scal_old,scal_new,aofs,dx,dt)
      implicit none
      include 'spec.h'
      integer nx

      real*8 scal_old(-1:nx  ,nscal)
      real*8 scal_new(-1:nx  ,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8 dx
      real*8 dt
      
      integer i,n

      n = Density
      do i = 0,nx-1
         scal_new(i,n) = scal_old(i,n) + dt * aofs(i,n)
      enddo
      
      scal_new(nx,n) = scal_new(nx-1,n)
      
      end
