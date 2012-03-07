      subroutine update_rho(scal_old,scal_new,aofs,dx,dt,time)
      implicit none
      include 'spec.h'
      real*8 scal_old(-2:nx+1,nscal)
      real*8 scal_new(-2:nx+1,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8 dx,dt,time

      integer i

      do i = 0,nx-1
         scal_new(i,Density) = scal_old(i,Density) + dt*aofs(i,Density)
      enddo
      
      call set_bc_s(scal_new,dx,time)

      end
