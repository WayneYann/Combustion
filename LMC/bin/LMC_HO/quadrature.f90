module quadrature

   implicit none
   
   public :: compute_integrals, compute_quadrature

contains

   subroutine compute_quadrature(q, f, m, dt)
      double precision, intent(out) :: q(:)
      double precision, intent(in ) :: f(nnodes,:)
      integer,          intent(in ) :: m
      double precision, intent(in ) :: dt
      
      if (m.eq.1) then
         q = (5.0*f(1,:)/24.0 + f(2,:)/3.0 - f(3,:)/24.0)*dt
      else
         q = (-f(1,:)/24.0 + f(2,:)/3.0 + 5*f(3,:)/24.0)*dt
      end if
      
   end subroutine compute_quadrature
   
   subroutine compute_integrals(I, a, d, r, dt, lo, hi)
      double precision, intent(out) :: I(nnodes-1,0:nx-1,nscal)
      double precision, intent(in ) :: a(nnodes,  0:nx-1,nscal)
      double precision, intent(in ) :: d(nnodes, -1:nx,  nscal)
      double precision, intent(in ) :: r(nnodes,  0:nx-1,nscal)
      double precision, intent(in ) :: dt
      
      double precision :: q1, q2
      double precision :: f1(nnodes,0:nx-1)
      integer :: n, m, i
      
      do n=1,nscal
         do i=lo,hi
            do m=1,nnodes
               f(m,i) = a(m,i,n) + d(m,i,n) + r(m,i,n)
            end do
         end do
         
         call compute_quadrature(I(1,:,n), f, 1, dt)
         call compute_quadrature(I(2,:,n), f, 2, dt)
      end do
      
   end subroutine compute_integrals

end module quadrature
