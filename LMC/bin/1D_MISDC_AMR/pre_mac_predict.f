      subroutine pre_mac_predict(vel_old,scal_old,gp,macvel,dx,dt,
     $                           lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8  vel_old(-2:nfine+1)
      real*8 scal_old(-2:nfine+1,nscal)
      real*8     gp(-1:nfine)
      real*8 macvel(0:nfine)
      real*8 dx, dt
      integer lo, hi, bc(2)
      
      real*8 slope(-1:nfine)
      real*8 dth
      real*8 dthx
      real*8 eps
      real*8 slo,shi
      integer i

      dth  = 0.5d0 * dt
      dthx = 0.5d0 * dt / dx
      eps = 1.d-6

      call mkslopes(vel_old,slope,lo,hi,bc)
      
      do i=lo,hi+1
         slo = vel_old(i-1) + (0.5d0 - dthx*vel_old(i-1))*slope(i-1) 
     $        - dth*gp(i-1)/scal_old(i-1,Density)
         shi = vel_old(i  ) - (0.5d0 + dthx*vel_old(i  ))*slope(i  )
     $        - dth*gp(i  )/scal_old(i  ,Density)
         if ( (slo+shi) .gt. eps) then
            macvel(i) = slo
         else if ( (slo+shi) .lt. -eps) then
            macvel(i) = shi
         else if ( (abs(slo+shi) .le. eps) .or. 
     $           (slo .le. 0.d0 .and. shi .ge. 0.d0)) then
            macvel(i) = 0.d0
         endif

         if (i .eq. lo .and. bc(1) .eq. 1) then
c     inflow
            macvel(i) = vel_old(i-1)
         end if

         if (i .eq. hi+1 .and. bc(2) .eq. 2) then
c     outflow
            macvel(i) = slo
         end if

      enddo

      end
