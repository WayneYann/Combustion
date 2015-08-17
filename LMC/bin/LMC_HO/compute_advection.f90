   subroutine get_diffdiff_terms(diffdiff, scal, beta, gamma_lo, gamma_hi, dx)

      implicit none
      include 'spec.h'

      double precision, intent(out) :: diffdiff( 0:nx-1)
      double precision, intent(in ) ::     scal(-2:nx+1,nscal)
      double precision, intent(in ) ::     beta(-2:nx+1,nscal)
      double precision, intent(in ) :: gamma_lo( 0:nx-1,Nspec)
      double precision, intent(in ) :: gamma_hi( 0:nx-1,Nspec)
      double precision, intent(in ) ::       dx

      real*8 dxsqinv,RWRK
      integer i,is,n,IWRK
      real*8 hm(Nspec,-1:nx)
      real*8 flux_lo(Nspec),flux_hi(Nspec)
      real*8 Y(Nspec,-1:nx)
      real*8 beta_lo, beta_hi, rho
   
      dxsqinv = 1.0/(dx*dx)

      diffdiff = 0

      do i=-1,nx
         rho = 0
!        compute density
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
!        compute Y = (rho*Y)/rho
         do n=1,Nspec
            Y(n,i) = scal(i,FirstSpec+n-1)/rho
         enddo
!        compute cell-centered h_m
         call CKHMS(scal(i,Temp),IWRK,RWRK,hm(1,i))
      end do

      do i=0,nx-1
         do n=1,Nspec
            is = FirstSpec + n - 1

!     compute -lambda/cp on faces
            if (coef_avg_harm.eq.1) then
               beta_lo = -2.d0 / (1.d0/beta(i,RhoH)+1.d0/beta(i-1,RhoH))
               beta_hi = -2.d0 / (1.d0/beta(i,RhoH)+1.d0/beta(i+1,RhoH))
            else
               beta_lo = -(beta(i  ,RhoH)+beta(i-1,RhoH)) /2.d0
               beta_hi = -(beta(i+1,RhoH)+beta(i  ,RhoH)) /2.d0
            end if

!     set face fluxes to -lambda/cp * grad Y_m
            flux_lo(n) = beta_lo*(Y(n  ,i) - Y(n,i-1))
            flux_hi(n) = beta_hi*(Y(n,i+1) - Y(n  ,i))

!     set face fluxes to h_m * (rho D_m - lambda/cp) grad Y_m
            flux_lo(n) = (flux_lo(n) + gamma_lo(i,n))*(hm(n,i-1)+hm(n,i))/2.d0
            flux_hi(n) = (flux_hi(n) + gamma_hi(i,n))*(hm(n,i+1)+hm(n,i))/2.d0
 
!     differential diffusion is divergence of face fluxes
            diffdiff(i) = diffdiff(i) + (flux_hi(n) - flux_lo(n))*dxsqinv

         end do
      end do
   end subroutine get_diffdiff_terms

   subroutine compute_advection(advection, scal_cc, vel, dx)
      use cell_conversions_module
      implicit none
      include 'spec.h'
      double precision, intent(out  ) :: advection( 0:nx-1,nscal)
      double precision, intent(in   ) ::   scal_cc(-2:nx+1,nscal)
      double precision, intent(in   ) ::       vel( 0:nx)
      double precision, intent(in   ) ::  dx
      
      double precision :: scal_face(0:nx, nscal)
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
            call cc_to_face(scal_face(:,n), scal_cc(:,n))
            !do i = 0,nx
            !   scal_face(i, n) = (-scal_cc(i-2,n) + 9.0*scal_cc(i-1,  n) &
            !                 + 9.0*scal_cc(i,  n) - scal_cc(i+1,n))/16.0
            !end do
            
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
             do i=0,nx-1
                advection(i,n) = - (vel(i+1)*scal_face(i+1,n) &
                                  - vel(i  )*scal_face(i  ,n)) / dx
             enddo
          endif
      enddo
   end subroutine compute_advection
