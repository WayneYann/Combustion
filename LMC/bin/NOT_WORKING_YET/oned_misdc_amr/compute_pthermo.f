      subroutine compute_pthermo(scal,lo,hi)
      implicit none
      include 'spec.h'
      real*8  scal(-2:nfine+1,nscal)
      integer lo,hi
      
      real*8 Y(Nspec), RWRK
      integer i,n
      integer ispec, IWRK

c     Define thermodynamic pressure
      do i=lo,hi
         do n = 1,nspec
            ispec = FirstSpec-1+n
            Y(n) = scal(i,ispec) / scal(i,Density)
         end do
         CALL CKPY(scal(i,Density),scal(i,Temp),Y,IWRK,RWRK,
     $             scal(i,RhoRT))
      end do
      scal(lo-1,RhoRT) = scal(lo,RhoRT)
      scal(hi+1,RhoRT) = scal(hi,RhoRT)
      
      end
