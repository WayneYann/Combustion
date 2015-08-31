   subroutine compute_pthermo(scal_cc)
      implicit none
      include 'spec.h'
      double precision, intent(inout) :: scal_cc(-2:nx+1,nscal)
      
      double precision :: Y_cc(Nspec)
      double precision :: rwrk
      integer          :: i, n, iwrk

      ! Define thermodynamic pressure
      do i=0,nx-1
         do n = 1,nspec
            Y_cc(n) = scal_cc(i,FirstSpec+n-1) / scal_cc(i,Density)
         end do
         
         CALL ckpy(scal_cc(i,Density),scal_cc(i,Temp),Y_cc,iwrk,rwrk,scal_cc(i,RhoRT))
      end do
   end