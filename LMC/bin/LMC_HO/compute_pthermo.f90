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
            Y(n) = scal(i,FirstSpec+n-1) / scal(i,Density)
         end do
         
         CALL CKPY(scal(i,Density),scal(i,Temp),Y,IWRK,RWRK,scal(i,RhoRT))
      end do
   end
