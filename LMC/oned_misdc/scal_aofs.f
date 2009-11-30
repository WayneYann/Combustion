      subroutine scal_aofs(scal_old,macvel,aofs,tforce,dx,dt)
      implicit none
      include 'spec.h'
      real*8 scal_old(-1:nx  ,nscal)
      real*8   macvel(0 :nx  )
      real*8   tforce(0 :nx-1,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8 dx
      real*8 time, dt
      
      real*8  sedge(0:nx  ,nscal)
      real*8  slope(0:nx-1)
      real*8 dth
      real*8 dthx
      real*8 eps
      real*8 slo,shi
      integer iconserv
      integer lo,hi,lo_lim,hi_lim
      integer ispec
      integer i,n
      integer use_bds

      real*8 Y(maxspec), RWRK, hmix
      integer IWRK
      
      dth  = 0.5d0 * dt
      dthx = 0.5d0 * dt / dx
      eps = 1.e-6
      
      use_bds = 0
      
      call set_bc_s(scal_old,dx,time)

      do n = 1,nscal
         if ((n.ne.Density) .and. (n.ne.RhoH) ) then
            
            if (n.eq.Temp) then
               iconserv = 0
            else
               iconserv = 1
            endif

            if (use_bds .eq. 1) then
               call bdsslopes(scal_old(-1,n),slope)
            else
               call  mkslopes(scal_old(-1,n),slope)
            endif
            do i = 1,nx-1
               slo = scal_old(i-1,n)+(0.5d0 - dthx*macvel(i))*slope(i-1)
               shi = scal_old(i  ,n)-(0.5d0 + dthx*macvel(i))*slope(i  )
               
               if (iconserv .eq. 1) then
                  slo = slo - dth * scal_old(i-1,n) * 
     $                 (macvel(i  )-macvel(i-1))/dx
                  shi = shi - dth * scal_old(i  ,n) * 
     $                 (macvel(i+1)-macvel(i  ))/dx
               endif
               
               slo = slo + dth*tforce(i-1,n)
               shi = shi + dth*tforce(i  ,n)
               
               if ( macvel(i) .gt. eps) then
                  sedge(i,n) = slo
               else if ( macvel(i) .lt. -eps) then
                  sedge(i,n) = shi
               else if ( abs(macvel(i)) .le. eps) then
                  sedge(i,n) = 0.5d0 * (slo + shi)
               endif
               
            enddo
            
            i = 0
            sedge(0,n) = scal_old(-1,n)
            
            i = nx
            sedge(i,n) = scal_old(i-1,n) + 
     $           (0.5d0 - dthx*macvel(i))*slope(i-1) 
            if (iconserv .eq. 1) then
               sedge(i,n) = sedge(i,n) - 
     $              dth * scal_old(i-1,n) * (macvel(i  )-macvel(i-1))/dx
            endif
            sedge(i,n) = sedge(i,n) + dth*tforce(i-1,n)
            
         endif
      enddo
      
c     Compute Rho on edges as sum of (Rho Y_i) on edges, rho.hmix as sum of (H_i.Rho.Y_i)
c     NOTE: Assumes Le=1 (no Le terms in RhoH equation)
      if (LeEQ1 .ne. 1) then
         print *,'Le != 1 terms not yet in aofs for RhoH'
      endif

      do i = 0,nx
         sedge(i,Density) = 0.d0
         do n = 1,nspec
            ispec = FirstSpec-1+n
            sedge(i,Density) = sedge(i,Density) + sedge(i,ispec)
         enddo
         do n = 1,nspec
            ispec = FirstSpec-1+n
            Y(n) = sedge(i,ispec) / sedge(i,Density)
         enddo
         call CKHBMS(sedge(i,Temp),Y,IWRK,RWRK,Hmix)
         sedge(i,RhoH) = Hmix * sedge(i,Density)
      enddo

      do n = 1,nscal
         if (n.eq.Temp) then
            do i = 0,nx-1
               aofs(i,n) = ( macvel(i+1)*sedge(i+1,n) 
     $                      -macvel(i  )*sedge(i  ,n)) / dx
               aofs(i,n) = aofs(i,n) - 
     $              (macvel(i+1)  - macvel(i)) * 0.5d0 * 
     $              ( sedge(i,n) + sedge(i+1,n)) /dx
             enddo
          else
             do i = 0,nx-1
                aofs(i,n) = ( macvel(i+1)*sedge(i+1,n) 
     $                       -macvel(i  )*sedge(i  ,n)) / dx
             enddo
          endif
          
c     Make these negative here so we can add as source terms later.
          do i = 0,nx-1
             aofs(i,n) = -aofs(i,n)
          enddo
       enddo
       
       end
