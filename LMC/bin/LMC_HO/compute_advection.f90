      subroutine compute_advection(advection, scal_cc, vel, dx, lo, hi)
      implicit none
      include 'spec.h'
      double precision, intent(out  ) :: advection(0 :nx-1,nscal)
      double precision, intent(in   ) ::   scal_cc(-2:nx+1,nscal)
      double precision, intent(in   ) ::       vel(0:nx-1)
      double precision, intent(in   ) ::  dx
      integer,          intent(in   ) ::  lo, hi
      
      
      double precision :: scal_face(0:nx, nscal)
      double precision :: RWRK
      integer          :: IWRK
      integer          :: i,n
      logical          :: compute_comp(nscal)
      
      ! we need to compute the advection term for: 
      ! density rho, density*enthalpy rho h, density*mass fractions rho Y_j
      ! (i.e. don't need to compute for RhoRT or temperature)
      do n = 1, nscal
         compute_comp(n) = .true.
      enddo
      compute_comp(RhoRT) = .false.
      compute_comp(Temp) = .false.
      
      ! we first compute face-values for the scalar quantities
      ! using the stencil, given cell-average values
      do n = 1,nscal
         if (compute_comp(n)) then
            do i = lo,hi+1
               scal_face(i, n) = (-scal_cc(i-2,n) + 9.0*scal_avg(i-1,  n) &
                             + 9.0*scal_cc(i,  n) - scal_avg(i+1,n))/16.0
            end do
            
!            call mkslopes(scal_old(:,n),slope,lo,hi,bc)
!            
!            do i=lo,hi+1
!               slo = scal_old(i-1,n)+(0.5d0)*slope(i-1)
!               shi = scal_old(i  ,n)-(0.5d0)*slope(i  )
!               
!               if ( macvel(i) .gt. eps) then
!                  sedge(i,n) = slo
!               else if ( macvel(i) .lt. -eps) then
!                  sedge(i,n) = shi
!               else if ( abs(macvel(i)) .le. eps) then
!                  sedge(i,n) = 0.5d0 * (slo + shi)
!               endif
!               
!               ! inflow
!               if (i .eq. lo) then
!                  sedge(i,n) = scal_old(i-1,n)
!               end if
!               
!               ! outflow
!               if (i .eq. hi+1) then
!                  sedge(i,n) = slo
!               end if
!               
!            enddo
            
         endif
      enddo
      
      ! compute the advection term
      do n = 1,nscal
         if (compute_comp(n)) then
             do i=lo,hi
                advection(i,n) = - (vel(i+1)*scal_face(i+1,n) &
                                  - vel(i  )*scal_face(i  ,n)) / dx
             enddo
          endif
      enddo
      
      end
