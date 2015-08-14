   
   ! use the mass fractions and temperature to compute the production rate wdot
   subroutine compute_production_rate(wdot, scal_cc)
      implicit none
      include 'spec.h'
      
      double precision, intent(out) ::    wdot( 0:nx-1, Nspec)
      double precision, intent(in ) :: scal_cc(-2:nx+1, nscal)
      
      double precision :: C(Nspec)
      double precision :: rwrk
      integer :: i, n, iwrk
      
      do i=0,nx-1
         do n=1,Nspec
            C(n) = scal_cc(i,FirstSpec+n-1)*invmwt(n)
         end do
         
         call CKWC(scal_cc(i,Temp), C, iwrk, rwrk, wdot(i,:))
      end do
   end subroutine compute_production_rate
