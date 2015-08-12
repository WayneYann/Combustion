      subroutine compute_pthermo(scal)
      implicit none
      include 'spec.h'
      double precision, intent(inout) :: scal(-2:nx+1,nscal)
      
      double precision :: Y(Nspec)
      double precision :: RWRK
      integer ::    i, n, IWRK

      ! Define thermodynamic pressure
      do i=0,nx-1
         do n = 1,nspec
            Y(n) = scal(i,FirstSpec-1+n) / scal(i,Density)
         end do
         
         CALL CKPY(scal(i,Density),scal(i,Temp),Y,IWRK,RWRK,scal(i,RhoRT))
      end do
      
      ! inflow
      scal(-1,RhoRT) = scal(0,RhoRT)

      ! outflow
      scal(nx,RhoRT) = scal(nx-1,RhoRT)
      
      end
