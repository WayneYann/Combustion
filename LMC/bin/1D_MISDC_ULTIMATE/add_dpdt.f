      subroutine add_dpdt(scal,pthermo,divu,umac,dx,dt,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8    scal(-2:nfine+1,nscal)
      real*8 pthermo(-2:nfine+1)
      real*8    divu( 0:nfine-1) ! divu normally has a ghost cell; be careful when passing in
      real*8    umac(0 :nfine  )
      real*8 dx
      real*8 dt
      integer lo,hi,bc(2)
      
      real*8 uadv,p_lo,p_hi
      real*8 ugradp
      real*8 dpdt
      integer i
      
      do i=lo,hi
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
 !        dpdt = dpdt - ugradp

         divu(i) = divu(i) + dpdt*dpdt_factor/pthermo(i)
      end do
      
      end


      subroutine add_dpdt_nodal(scal,pthermo,divu,vel,dx,dt,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8    scal(-2:nfine+1,nscal)
      real*8 pthermo(-2:nfine+1)
      real*8    divu(0:nfine-1)
      real*8     vel(-2:nfine+1)
      real*8 dx
      real*8 dt
      integer lo,hi,bc(2)
      
      real*8 uadv,p_lo,p_hi
      real*8 ugradp
      real*8 dpdt
      integer i
      
      do i=lo,hi
         uadv = vel(i)

         if (i .eq. lo .and. bc(1) .eq. 1) then
!     inflow
            p_lo = pthermo(i)
            p_hi = pthermo(i+1)
         else if (i .eq. hi .and. bc(2) .eq. 2) then
!     outflow
            p_lo = pthermo(i-1)
            p_hi = pthermo(i)
         else
            if (uadv .ge. 0) then
               p_lo = pthermo(i-1)
               p_hi = pthermo(i)
            else                  
               p_lo = pthermo(i)
               p_hi = pthermo(i+1)
            end if
         end if

         ugradp = uadv * (p_hi - p_lo) / dx 
         dpdt = (pthermo(i) - Pcgs) / dt
!         dpdt = dpdt - ugradp
         divu(i) = divu(i) + dpdt*dpdt_factor/pthermo(i)
      end do
      
      end
