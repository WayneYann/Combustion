      subroutine bdsslopes(nx,s,slope)
      implicit none
      
      integer nx
      real*8     s(-1:nx  )
      real*8 slope( 0:nx-1)
      
      real*8 sint(0:nx)
      real*8 slo,shi
      real*8 smin,smax
      integer i

      real*8 sumloc,redfac
      real*8 diff_lo,diff_hi
      real*8 sumdif,sgndif
      real*8 smax_lo,smax_hi
      real*8 smin_lo,smin_hi
      real*8 sc_lo,sc_hi
      real*8 eps
      integer kdp
      integer ll
      integer unlim

      eps = 1.d-8

      unlim = 0

      do i = 1,nx-1
         sint(i) =       ( -(s(i-2) + s(i+1)) 
     $        + 7.d0 * (s(i-1) + s(i  )) ) / 12.d0
      enddo

      sint( 0) = 0.5d0 * (s(0 ) + s(  -1))
      sint(nx) = 0.5d0 * (s(nx) + s(nx-1))

      do i = 0,nx-1
         slope(i) = sint(i+1) - sint(i)
      enddo

      if (unlim .eq. 0) then
c     Now apply min/max limiting.
         do i = 0,nx-1
            shi = s(i) + 0.5d0*slope(i)
            smax = max(s(i),s(i+1))
            smin = min(s(i),s(i+1))
            if (shi .gt. smax)
     $           slope(i) = 2.d0*(smax-s(i))
            if (shi .lt. smin) 
     $           slope(i) = 2.d0*(smin-s(i))
            
            slo = s(i) - 0.5d0*slope(i)
            smax = max(s(i),s(i-1))
            smin = min(s(i),s(i-1))
            if (slo .gt. smax)
     $           slope(i) = 2.d0*(s(i)-smax)
            if (slo .lt. smin) 
     $           slope(i) = 2.d0*(s(i)-smin)

         enddo
      endif

c       Traditional limiting (gives same result)
c       do i = 0,nx-1

c         smin_hi = min(s(i),s(i+1))
c         smin_lo = min(s(i),s(i-1))
c         smax_hi = max(s(i),s(i+1))
c         smax_lo = max(s(i),s(i-1))
c
c         sc_hi = s(i) + 0.5d0 * slope(i)
c         sc_lo = s(i) - 0.5d0 * slope(i)
c
c         sc_lo = max(min(sc_lo, smax_lo), smin_lo)
c         sc_hi = max(min(sc_hi, smax_hi), smin_hi)
c
c         do ll = 1,1 
c
c           sumloc = 0.5d0*(sc_lo + sc_hi)
c           sumdif = (sumloc - s(i))*2.d0
c           sgndif = sign(1.d0,sumdif)
c
c           diff_lo = (sc_lo - s(i))*sgndif
c           diff_hi = (sc_hi - s(i))*sgndif
c
c           kdp = 0
c           if (diff_lo.gt.eps) kdp = kdp+1
c           if (diff_hi.gt.eps) kdp = kdp+1
c
c           First adjust lo value if it violates min or max.
c           if (diff_lo .gt. eps) then
c             redfac = sumdif * sgndif / dble(kdp)
c             if (sgndif .gt. 0.d0) then
c                redfac = min(redfac,(sc_lo-smin_lo))
c             else
c                redfac = min(redfac,(smax_lo-sc_lo))
c             endif
c             sumdif = sumdif - redfac*sgndif
c             sc_lo = sc_lo - redfac*sgndif
c             kdp = kdp - 1
c           endif
c
c           Then adjust hi value
c           if (diff_hi .gt. eps) then
c             redfac = sumdif * sgndif / dble(kdp)
c             if (sgndif .gt. 0.d0) then
c                redfac = min(redfac,(sc_hi-smin_hi))
c             else
c                redfac = min(redfac,(smax_hi-sc_hi))
c             endif
c             sumdif = sumdif - redfac*sgndif
c             sc_hi = sc_hi - redfac*sgndif
c             kdp = kdp - 1
c           endif

c           slope(i) = sc_hi - sc_lo

c         enddo
c       enddo

      end
