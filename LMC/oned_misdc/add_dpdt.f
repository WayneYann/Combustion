      subroutine add_dpdt(pthermo,divu,umac,dx,dt)
      implicit none
      include 'spec.h'
      real*8 pthermo(-1:nx)
      real*8  divu(0 :nx-1)
      real*8  umac(0 :nx  )
      real*8 dx
      real*8 dt
      
      real*8 uadv,p_lo,p_hi
      real*8 ugradp
      real*8 denom
      real*8 dpdt
      real*8 dpdt_max
      integer i,n
      integer ispec
      
      dpdt_max = 0.d0
      do i = 0,nx-1
         uadv = 0.5d0 * (umac(i) + umac(i+1))
         if (umac(i) .ge. 0.d0) then
            p_lo = pthermo(i-1)
         else
            p_lo = pthermo(i)
         endif
         if (umac(i+1) .ge. 0.d0) then
            p_hi = pthermo(i)
         else
            p_hi = pthermo(i+1)
         endif
         ugradp = uadv * (p_hi - p_lo) / dx 
         dpdt = (pthermo(i) - Pcgs) / dt
         dpdt = dpdt - ugradp        
c     FIXME: compute gamma here
         dpdt = dpdt / 1.4d0
         denom = MIN(pthermo(i),Pcgs)
         dpdt = dpdt/denom
         dpdt = dpdt_factor * dpdt
         divu(i) = divu(i) + dpdt
         dpdt_max = MAX(dpdt_max,ABS(dpdt))
      end do
      
      write(6,1000) dpdt_max
 1000 format('DPDT norm     = ',f14.4)
      
      end
