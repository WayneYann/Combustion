      subroutine update_rho(nx,scal_old,scal_new,aofs,dx,dt)
      implicit none
      include 'spec.h'
      integer nx

      real*8 scal_old(-1:nx  ,nscal)
      real*8 scal_new(-1:nx  ,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8 dx,dt

      integer i

      do i = 0,nx-1
         scal_new(i,Density) = scal_old(i,Density) + dt*aofs(i,Density)
      enddo
      end
