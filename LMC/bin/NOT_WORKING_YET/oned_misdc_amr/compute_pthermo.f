      subroutine compute_pthermo(scal,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8  scal(-2:nfine+1,nscal)
      integer lo,hi
      integer bc(2)
      
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

      if (bc(1) .eq. 1) then
c     inflow
         scal(lo-1,RhoRT) = scal(lo,RhoRT)
      end if

      if (bc(2) .eq. 2) then
c     outflow
         scal(hi+1,RhoRT) = scal(hi,RhoRT)
      end if
      
      end
