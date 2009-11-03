      subroutine mkslopes(nx,s,slope)
      implicit none
      integer nx
      integer ncomp
      real*8     s(-1:nx  )
      real*8 slope( 0:nx-1)

      real*8 slxscr(-1:nx,4)
      real*8 slo,shi,slim,smid,ds
      real*8 smax,smin
      integer i
      integer slope_order
      integer unlim

      integer cen,lim,flag,fromm
      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flag = 3 )
      parameter( fromm = 4 )

      slope_order = 4
      unlim       = 0
      
      if (slope_order .eq. 1) then
         do i = 0,nx-1
            slope(i) = 0.d0
         enddo

c     2nd order slopes
      elseif (slope_order .eq. 2) then
         
         if (unlim .eq. 0) then
            do i = 0,nx-1
               shi = 2.0d0*(s(i+1) - s(i  ))
               slo = 2.0d0*(s(i  ) - s(i-1))
               smid = 0.5d0*(s(i+1) - s(i-1))
               if (i.eq.0) smid = (s(1)+3.d0*s(0)-4.d0*s(-1))/3.d0
               slim = min(abs(slo),abs(shi))
               slim = min(abs(smid),slim)
               if (slo*shi .lt. 0.d0) then
                  slope(i) = 0.d0
               else
                  slope(i) = sign(1.d0,smid)*slim
               endif
            enddo
            slope(nx-1) = 0.d0

         else
c     Unlimited slopes
            do i = 0,nx-1
               slope(i) = 0.5d0*(s(i+1) - s(i-1))
            enddo
            slope(nx-1) = 0.d0
         endif

c     4th order slopes
      else
         
         if (unlim .eq. 0) then
            do i = 1,nx-1
               shi = 2.0d0*(s(i+1) - s(i  ))
               slo = 2.0d0*(s(i  ) - s(i-1))
               slxscr(i,cen) = 0.5d0*(s(i+1)-s(i-1))
               
               slxscr(i,lim)  = min(abs(slo),abs(shi))
               if (slo*shi .le. 0.d0) then
                  slxscr(i,lim)  = 0.d0
               endif
               slxscr(i,flag) = sign(1.d0,slxscr(i,cen))
               slxscr(i,fromm)= slxscr(i,flag)*
     &              min(slxscr(i,lim),abs(slxscr(i,cen)))

            enddo
c     FIXME
c     Inflow (EXT_DIR) boundary value specified at i=0
            i = 0
            smid = -16.d0/15.d0*s(-1) + 0.5d0*s(0) 
     &           + 2.d0/3.d0*s(1) - 0.1d0*s(2)
            shi = 2.0d0*(s(1) - s( 0))
            slo = 2.0d0*(s(0) - s(-1))
            slim = min(abs(slo),abs(shi))
            slim = min(abs(smid),slim)
            if (slo*shi .lt. 0.d0) then
               slope(0) = 0.d0
            else
               slope(0) = sign(1.d0,smid)*slim
            endif

            slxscr(0,flag) = sign(1.d0,slxscr(0,cen))
            slxscr(0,fromm) = slope(0)
            
            do i = 1,nx-1
               ds = (4.d0/3.d0) * slxscr(i,cen) -
     &              (1.d0/6.d0) * (slxscr(i+1,fromm)+slxscr(i-1,fromm))
               slope(i) = slxscr(i,flag)*min(abs(ds),slxscr(i,lim))
            enddo

         else
c     Unlimited slopes
            do i = 0,nx-1
               slxscr(i,cen) = 0.5d0*(s(i+1)-s(i-1))
            enddo
            do i = 1,nx-2
               slope(i) = (4.d0/3.d0)*slxscr(i,cen)
     &              -(1.d0/6.d0)*(slxscr(i+1,cen) + slxscr(i-1,cen))
            enddo
            slope(0   ) = 0.d0
         endif

         slope(nx-1) = 0.d0

      endif

      end
