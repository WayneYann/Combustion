      subroutine compute_production_rate(wdot, scal_cc, lo, hi)
      implicit none
      include 'spec.h'
      
      double precision, intent(out) ::    wdot( 0:nx-1, Nspec)
      double precision, intent(in ) :: scal_cc(-2:nx+1, nscal)
      integer,          intent(in ) ::  lo, hi
      
      double precision :: C(Nspec)
      double precision :: rwrk
      integer :: i, n, iwrk
      
      do i=lo,hi
         do n=1,Nspec
            C(n) = scal_cc(i,NSpec+n-1)*invmwt(n)
         end do
         
         call CKWC(scal_cc(i,Temp), C, iwrk, rwrk, wdot(i,:))
      end do
      
      end subroutine compute_production_rate
