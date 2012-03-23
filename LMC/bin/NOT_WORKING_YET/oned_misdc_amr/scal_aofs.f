      subroutine scal_aofs(scal_old,macvel,aofs,divu,tforce,
     $                     dx,dt,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8 scal_old(-2:nfine+1,nscal)
      real*8   macvel(0 :nfine  )
      real*8   tforce(-1:nfine  ,nscal)
      real*8     aofs(0 :nfine-1,nscal)
      real*8     divu(-1:nfine)
      real*8 dx
      real*8 dt
      integer lo,hi,bc(2)
      
      real*8  sedge(0 :nfine,nscal)
      real*8  slope(-1:nfine)
      real*8 dth
      real*8 dthx
      real*8 eps
      real*8 slo,shi
      integer iconserv
      integer ispec
      integer i,n

      real*8 Y(Nspec), RWRK, hmix
      integer IWRK

      logical compute_comp(nscal)
      
      dth  = 0.5d0 * dt
      dthx = 0.5d0 * dt / dx
      eps = 1.e-6
      
      do n = 1, nscal
         compute_comp(n) = .true.
      enddo
      compute_comp(Density) = .false.

      if (use_strang) then
         ! predict everything
      else
         ! no need to predict temperature
         compute_comp(Temp) = .false.
      end if

      do n = 1,nscal
         if (compute_comp(n) ) then
            
            if (n.eq.Temp) then
               iconserv = 0
            else
               iconserv = 1
            endif

            call mkslopes(scal_old(:,n),slope,lo,hi,bc)

            do i=lo+1,hi
               slo = scal_old(i-1,n)+(0.5d0 - dthx*macvel(i))*slope(i-1)
               shi = scal_old(i  ,n)-(0.5d0 + dthx*macvel(i))*slope(i  )
               
               if (iconserv .eq. 1) then
                  slo = slo - dth * scal_old(i-1,n) * divu(i-1)
                  shi = shi - dth * scal_old(i  ,n) * divu(i)
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
            
            i = lo
            sedge(0,n) = scal_old(lo-1,n)
            
            i = hi+1
            sedge(i,n) = scal_old(i-1,n) + 
     $           (0.5d0 - dthx*macvel(i))*slope(i-1) 
            if (iconserv .eq. 1) then
               sedge(i,n) = sedge(i,n) - 
     $              dth * scal_old(i-1,n) * divu(i-1)
            endif
            sedge(i,n) = sedge(i,n) + dth*tforce(i-1,n)
            
         endif
      enddo
      
      do i=lo,hi+1
         sedge(i,Density) = 0.d0
c        compute Rho on edges as sum of (Rho Y_i) on edges, 
         do n = 1,nspec
            ispec = FirstSpec-1+n
            sedge(i,Density) = sedge(i,Density) + sedge(i,ispec)
         enddo
c        compute rho.hmix as sum of (H_i.Rho.Y_i)
         if ( .not.compute_comp (RhoH) ) then 
            do n = 1,nspec
               ispec = FirstSpec-1+n
               Y(n) = sedge(i,ispec) / sedge(i,Density)
            enddo
            call CKHBMS(sedge(i,Temp),Y,IWRK,RWRK,Hmix)
            sedge(i,RhoH) = Hmix * sedge(i,Density)
         endif
      enddo

      do n = 1,nscal
         if (n.eq.Temp) then
            do i=lo,hi
               aofs(i,n) = ( macvel(i+1)*sedge(i+1,n) 
     $                      -macvel(i  )*sedge(i  ,n)) / dx
               aofs(i,n) = aofs(i,n) - 
     $              (macvel(i+1)  - macvel(i)) * 0.5d0 * 
     $              ( sedge(i,n) + sedge(i+1,n)) /dx
             enddo
          else
             do i=lo,hi
                aofs(i,n) = ( macvel(i+1)*sedge(i+1,n) 
     $                       -macvel(i  )*sedge(i  ,n)) / dx
             enddo
          endif
          
c     Make these negative here so we can add as source terms later.
          do i=lo,hi
             aofs(i,n) = -aofs(i,n)
          enddo
       enddo
       
       end

