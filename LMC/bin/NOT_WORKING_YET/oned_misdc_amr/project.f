      subroutine project_level(vel_new,rhohalf,divu,
     $                         press_old,press_new,dx,dt,lo,hi,bc)
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
         press_new(lo-1) = press_new(lo)
         press_new(hi+2) = press_new(hi+1)

      endif

c     also need to set outflow boundary
      call set_bc_v(vel_new,lo,hi,bc)

      end


      subroutine project_composite(vel,divu,dx,dt,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8   vel(0:nlevs-1,-2:nfine+1)
      real*8  divu(0:nlevs-1, 0:nfine)
      real*8    dx(0:nlevs-1)
      real*8    dt(0:nlevs-1)
      integer   lo(0:nlevs-1),hi(0:nlevs-1)
      integer   bc(0:nlevs-1,2)

c     local
      integer i
      real*8 divu_node,vel_temp,offset

      if (nlevs .gt. 2) then
         print*,"ERROR: composite_project assumes nlevs <= 2"
         stop
      end if

      call set_bc_v(vel(0,:),lo(0),hi(0),bc(0,:))

c     integrate the coarse level first as if fine level didn't exist
      vel(0,0) = vel(0,-1) + 0.5*divu(0,0)*dx(0)
      do i=1,hi(0)
         divu_node = (divu(0,i) + divu(0,i-1))*0.5d0
         vel(0,i) = vel(0,i-1) + divu_node*dx(0)
      end do

      if (nlevs .eq. 2) then
c     integrate the fine level
         vel(1,lo(1)) = vel(0,lo(1)/2-1) 
     $        + 0.5d0*dx(0)*divu(0,lo(1)/2-1)
     $        + 0.5d0*dx(1)*divu(1,lo(1))
         do i=lo(1)+1,hi(1)
            divu_node = (divu(1,i) + divu(1,i-1))*0.5d0
            vel(1,i) = vel(1,i-1) + divu_node*dx(1)
         end do
         
c     compute mismatch by first computing what the next coarse
c     cell would have been
         vel_temp = vel(1,hi(1))
     $        + 0.5d0*dx(1)*divu(1,hi(1))
     $        + 0.5d0*dx(0)*divu(0,(hi(1)+1)/2)

         offset = vel_temp - vel(0,(hi(1)+1)/2)

c     then offset the remainder of the coarse solution
         do i=(hi(1)+1)/2,hi(0)
            vel(0,i) = vel(0,i) + offset
         end do

      end if

c     fill outflow ghost cells
      call set_bc_v(vel(0,:),lo(0),hi(0),bc(0,:))

c     fillpatch
      if (nlevs .eq. 2) then
         call fill_ghost(vel,vel,1,lo,hi,0.d0,0.d0,dx,dt)
      end if

      end
