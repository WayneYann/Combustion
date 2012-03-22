      subroutine macproj(macvel,divu,dx,lo,hi)
      implicit none
      include 'spec.h'
      real*8 macvel(0 :nfine)
      real*8   divu(0 :nfine-1)
      real*8 dx
      integer lo,hi

      integer i
      print *,'... mac_projection'

c     macvel(lo) was set in pre_mac_predict and remains unchanged
      
      do i=lo+1,hi+1
         macvel(i) = macvel(i-1) + divu(i-1)*dx
      end do
      
      end
