      subroutine mkslopes(scal,slope)
      implicit none
      include 'spec.h'
      real*8  scal(-2:nx+1)
      real*8 slope( 0:nx-1)

      real*8 slo,shi,slim,smid
      integer i

      integer cen,lim,flag,fromm
      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flag = 3 )
      parameter( fromm = 4 )

      if (unlim .eq. 0) then
         do i = 0,nx-1
            shi = 2.0d0*(scal(i+1) - scal(i  ))
            slo = 2.0d0*(scal(i  ) - scal(i-1))
            smid = 0.5d0*(scal(i+1) - scal(i-1))
            if (i.eq.0) smid = (scal(1)+3.d0*scal(0)-4.d0*scal(-1))/3.d0
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
         do i = 0,nx-1
            slope(i) = 0.5d0*(scal(i+1) - scal(i-1))
         enddo
         slope(nx-1) = 0.d0
      endif

      end
