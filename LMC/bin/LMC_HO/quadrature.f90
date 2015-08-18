module quadrature_module

   use cell_conversions_module
   implicit none
   private
   include 'spec.h'
   public :: compute_integrals_avg, compute_integrals_cc, compute_quadrature

contains
   
   ! return the integral of f over a subinterval indicated by m
   ! m = 1, return the integral over the interval [0,  dt/2]
   ! m = 2, return the integral over the interval [dt/2, dt]
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
   
   ! compute the cell-average value of the integrals from 
   ! 0 to dt/2 and from dt/2 to dt of A+D+R
   ! A and D are passed in as cell-averaged quantities
   ! and R is passed in as a cell-centered quantity
   subroutine compute_integrals_avg(q, a, d, r, dt)
      double precision, intent(out) :: q(nnodes-1,0:nx-1,nscal)
      double precision, intent(in ) :: a(nnodes,  0:nx-1,nscal)
      double precision, intent(in ) :: d(nnodes,  0:nx-1,nscal)
      double precision, intent(in ) :: r(nnodes,  0:nx-1,Nspec)
      double precision, intent(in ) :: dt
      
      double precision :: r_avg(nnodes, 0:nx-1, nscal)
      double precision :: f(nnodes,0:nx-1)
      integer :: n, m, i
      
      ! convert the R term to cell-averaged
      r_avg = 0
      do n=1,Nspec
         do m=1,nnodes
            call extrapolate_cc_to_avg(r_avg(m,:,FirstSpec+n-1), r(m,:,n))
         end do
      end do
      
      ! compute the quadratures
      do n=1,nscal
         do i=0,nx-1
            do m=1,nnodes
               f(m,i) = a(m,i,n) + d(m,i,n) + r_avg(m,i,n)
            end do
            call compute_quadrature(q(1,i,n), f(:,i), 1, dt)
            call compute_quadrature(q(2,i,n), f(:,i), 2, dt)
         end do
      end do
      
   end subroutine compute_integrals_avg
   
   ! compute the cell-centered value of the integrals from 
   ! 0 to dt/2 and from dt/2 to dt of A+D+R
   ! A and D are passed in as cell-averaged quantities
   ! and R is passed in as a cell-centered quantity
   subroutine compute_integrals_cc(q, a, d, r, dt)
      double precision, intent(out) :: q(nnodes-1,0:nx-1,nscal)
      double precision, intent(in ) :: a(nnodes,  0:nx-1,nscal)
      double precision, intent(in ) :: d(nnodes,  0:nx-1,nscal)
      double precision, intent(in ) :: r(nnodes,  0:nx-1,Nspec)
      double precision, intent(in ) :: dt
      
      double precision :: a_cc(nnodes, 0:nx-1,nscal)
      double precision :: d_cc(nnodes, 0:nx-1, nscal)
      double precision :: r_cc(nnodes, 0:nx-1, nscal)
      
      double precision :: f(nnodes,0:nx-1)
      integer :: n, m, i
      
      ! convert A and D to cell-centered quantities
      do n=1,nscal
         do m=1,nnodes
            call extrapolate_avg_to_cc(a_cc(m,:,n), a(m,:,n))
            call extrapolate_avg_to_cc(d_cc(m,:,n), d(m,:,n))
         end do
      end do
      
      ! R is zero for all scalars except mass fractions
      r_cc = 0
      do n=1,Nspec
         do m=1,nnodes
            r_cc(m,:,FirstSpec+n-1) = r(m,:,n)
         end do
      end do
      
      ! compute the quadratures
      do n=1,nscal
         do i=0,nx-1
            do m=1,nnodes
               f(m,i) = a_cc(m,i,n) + d_cc(m,i,n) + r_cc(m,i,n)
            end do
            call compute_quadrature(q(1,i,n), f(:,i), 1, dt)
            call compute_quadrature(q(2,i,n), f(:,i), 2, dt)
         end do
      end do
      
   end subroutine compute_integrals_cc

end module quadrature_module
