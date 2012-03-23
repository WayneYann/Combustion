      subroutine project(vel_new,rhohalf,divu,
     $                   press_old,press_new,dx,dt,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8   vel_new(-2:nfine+1)
      real*8   rhohalf(0 :nfine-1)
      real*8      divu(-1:nfine)
      real*8 press_old(-1:nfine+1)
      real*8 press_new(-1:nfine+1)
      real*8 dx
      real*8 dt
      integer lo,hi,bc(2)
      
      integer i
      real*8 phi(0:nfine)
      real*8 vel_star(0:nfine-1)
      real*8 divu_node
      real*8 gp

      ! if dt < 0, we simply project vel_new and leave pressure
      ! to zero.  This should only be the case in the initial
      ! projection and the divu_iters.

      if (dt .gt. 0) then         

         ! project v_star = vel_new + dt*gp^old/rho
         do i=lo,hi
            gp = (press_old(i+1)-press_old(i))/dx
            vel_star(i) = vel_new(i) + dt * gp / rhohalf(i)
         end do

      endif

c     Build vel_new directly, since we have bc and div(v)=s
c     Get boundary value, vel(-1) at inlet wall, and integrate explicitly.
      call set_bc_v(vel_new,lo,hi,bc)
      vel_new(0) = vel_new(-1) + 0.5*divu(0)*dx
      do i=lo+1,hi
         divu_node = (divu(i) + divu(i-1))*0.5d0
         vel_new(i) = vel_new(i-1) + divu_node*dx
      end do
      
      if (dt .gt. 0.) then

c     For this projection, vel_new = vel_star - dt*gp^new/rho
c     therefore gp = -(vel_new-vel_star)/dt * rhohalf
         phi(hi+1) = 0.d0
         do i=hi,lo,-1
            gp = -(vel_new(i) - vel_star(i))/dt * rhohalf(i)
            phi(i) = phi(i+1) - gp*dx
         enddo
         do i=lo,hi+1
            press_new(i) = phi(i)
         enddo

      endif

c     also need to set outflow boundary
      call set_bc_v(vel_new,lo,hi,bc)

      end
