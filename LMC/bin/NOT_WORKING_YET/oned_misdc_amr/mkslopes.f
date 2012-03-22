      subroutine mkslopes(scal,slope,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8  scal(-2:nfine+1)
      real*8 slope(-1:nfine  )
      integer lo,hi,bc(2)

      real*8 slo,shi,slim,smid
      integer i

      integer cen,lim,flag,fromm
      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flag = 3 )
      parameter( fromm = 4 )

      if (unlim .eq. 0) then
         do i=lo-1,hi+1
            shi = 2.0d0*(scal(i+1) - scal(i  ))
            slo = 2.0d0*(scal(i  ) - scal(i-1))
            smid = 0.5d0*(scal(i+1) - scal(i-1))
            if (i .eq. lo .and. bc(1) .eq. 1) then
               ! inflow
               ! value in ghost cell is value at inflow face
               smid = (scal(1)+3.d0*scal(0)-4.d0*scal(-1))/3.d0
            end if
            slim = min(abs(slo),abs(shi))
            slim = min(abs(smid),slim)
            if (slo*shi .lt. 0.d0) then
               slope(i) = 0.d0
            else
               slope(i) = sign(1.d0,smid)*slim
            endif
         enddo
      else
         do i=lo-1,hi+1
            slope(i) = 0.5d0*(scal(i+1) - scal(i-1))
         enddo
      endif

      if (bc(1) .eq. 1) then
         ! inflow
         slope(lo-1) = 0.d0
      end if

      if (bc(2) .eq. 2) then
         ! outflow
         slope(hi  ) = 0.d0
         slope(hi+1) = 0.d0
      end if

      end
