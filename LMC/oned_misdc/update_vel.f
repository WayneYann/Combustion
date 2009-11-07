      subroutine update_vel(nx,vel_old,vel_new,gp,rhohalf,
     $                      macvel,sedge,dx,dt)
      implicit none
c      include 'spec.h'

      integer nx

      real*8  vel_old(-1:nx)
      real*8  vel_new(-1:nx)
      real*8       gp(0 :nx-1)
      real*8  rhohalf(0 :nx-1)
      real*8   macvel(0 :nx)
      real*8    sedge(0 :nx)
      real*8 dx
      real*8 dt
      
      real*8 aofs
      
      integer i,n
      
      do i = 0,nx-1
         
         aofs = ( (macvel(i+1)*sedge(i+1) - macvel(i)*sedge(i))  -
     $        (macvel(i+1)-macvel(i)) * 
     $        0.5d0 * (sedge(i)+sedge(i+1)) ) /dx
         
         vel_new(i) = vel_old(i) + dt * ( aofs - gp(i) / rhohalf(i) )
      enddo
      vel_new(nx) = vel_new(nx-1)
      
      end
