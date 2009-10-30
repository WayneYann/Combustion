        subroutine vel_edge_states(nx,vel_old,scal_old,gp,
     $                             macvel,sedge,dx,dt)

        implicit none

        include 'nums.fi'
        include 'sndata.fi'

c       Quantities passed in 
        integer nx
        real*8  vel_old(-1:nx)
        real*8 scal_old(-1:nx,nscal)
        real*8       gp(0:nx-1)
        real*8   macvel(0:nx)
        real*8   sedge(0:nx)
        real*8 dx
        real*8 dt

c       Local variables
        real*8 slope(0:nx-1)
        real*8 dth
        real*8 dthx
        real*8 eps
        real*8 slo,shi
        integer i,n

        dth  = 0.5d0 * dt
        dthx = 0.5d0 * dt / dx
        eps = 1.e-6

        call mkslopes(nx,vel_old,slope)

        do i = 1,nx-1
          slo = vel_old(i-1) + (0.5 - dthx*vel_old(i-1))*slope(i-1)
     $        - dth*gp(i-1)/scal_old(i-1,Density)
          shi = vel_old(i  ) - (0.5 + dthx*vel_old(i))*slope(i  )
     $        - dth*gp(i  )/scal_old(i  ,Density)
          if ( macvel(i) .gt. eps) then
            sedge(i) = slo
          else if ( macvel(i) .lt. -eps) then
            sedge(i) = shi
          else if ( abs(macvel(i)) .le. eps) then
            sedge(i) = 0.5d0 * (slo + shi)
          endif
         enddo

         i = 0
         sedge(0) = vel_old(-1)

         i = nx
         sedge(i) = vel_old(i-1) + 
     $      (0.5d0 - dthx*vel_old(i-1))*slope(i-1) 
     $     - dth*gp(i-1)/scal_old(i-1,Density)

        end
