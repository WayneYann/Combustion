      subroutine compute_pthermo(scal, lo, hi)
      implicit none
      include 'spec.h'
      double precision, intent(inout) :: scal(-2:nx+1,nscal)
      integer,          intent(in   ) :: lo, hi
      
      double precision :: Y(Nspec)
      double precision :: RWRK
      integer ::    i, n, IWRK

      ! Define thermodynamic pressure
      do i=lo,hi
         do n = 1,nspec
            Y(n) = scal(i,FirstSpec-1+n) / scal(i,Density)
         end do
         
         CALL CKPY(scal(i,Density),scal(i,Temp),Y,IWRK,RWRK,scal(i,RhoRT))
      end do
      
      ! inflow
      scal(lo-1,RhoRT) = scal(lo,RhoRT)

      ! outflow
      scal(hi+1,RhoRT) = scal(hi,RhoRT)
      
      end
