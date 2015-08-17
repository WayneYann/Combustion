module quadrature_module

   use cell_conversions_module
   implicit none
   private
   include 'spec.h'
   public :: compute_integrals, compute_quadrature

contains
   
   subroutine compute_quadrature(q, f, m, dt)
      double precision, intent(out) :: q
      double precision, intent(in ) :: f(nnodes)
      integer,          intent(in ) :: m
      double precision, intent(in ) :: dt
      
      if (m.eq.1) then
         q = (5*f(1) + 8*f(2) - f(3))*dt/24.0
      else
         q = (-f(1) + 8*f(2) + 5*f(3))*dt/24.0
      end if
      
   end subroutine compute_quadrature
   
   subroutine compute_integrals(q, a, d, r, dt)
      double precision, intent(out) :: q(nnodes-1,0:nx-1,nscal)
      double precision, intent(in ) :: a(nnodes,  0:nx-1,nscal)
      double precision, intent(in ) :: d(nnodes,  0:nx-1,nscal)
      double precision, intent(in ) :: r(nnodes,  0:nx-1,Nspec)
      double precision, intent(in ) :: dt
      
      double precision :: r_avg(nnodes, 0:nx-1, nscal)
      double precision :: f(nnodes,0:nx-1)
      integer :: n, m, i
      
      r_avg = 0
      do n=1,Nspec
         do m=1,nnodes
            call extrapolate_cc_to_avg(r_avg(m,:,FirstSpec+n-1), r(m,:,n))
         end do
      end do
      
      do n=1,nscal
         do i=0,nx-1
            do m=1,nnodes
               f(m,i) = a(m,i,n) + d(m,i,n) + r_avg(m,i,n)
            end do
            call compute_quadrature(q(1,i,n), f(:,i), 1, dt)
            call compute_quadrature(q(2,i,n), f(:,i), 2, dt)
         end do
      end do
      
   end subroutine compute_integrals

end module quadrature_module
