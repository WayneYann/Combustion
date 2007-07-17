c#undef BL_LANG_CC
c#ifndef BL_LANG_FORT
c#define BL_LANG_FORT
c#endif

c#define BL_SPACEDIM 2

c#include "ArrayLim.H"

c     evaluate the "operator" (update and fluxes) for the Runge-Kutta 
c     implementation of the diffusion step using multicomponent coefficients
c     with Soret and Dufour effects

c     to do:
c     1) must adjust diffusion coefficients by thickness factors
c     2) unsure of LMC support for r-z coordinates

c     arguments are alphabetical, mostly:
c
c     domain_lo, domain_hi,             ! INPUT limits of valid region of the domain
c     lo, hi,                           ! INPUT limits of valid region of the box
c     areax, DIMS(areax),               ! INPUT areas of the faces perendicular to x axis
c     areay, DIMS(areay),               ! INPUT areas of the faces perpendicular to y axis
c     bc,                               ! INPUT boundary condition array for all comps
c     dt,                               ! INPUT timestep
c     dx,                               ! INPUT physical dimensions of grid cells
c     index_of_firstY,                  ! INPUT index of rho Y for the first species in the state
c     index_of_lastY,                   ! INPUT index of rho Y for the last species in the state
c     index_of_rho,                     ! INPUT index of rho in the state
c     index_of_rhoH,                    ! INPUT index of rho H in the state
c     index_of_T,                       ! INPUT index of T in the state
c     maximum_error,                    ! INPUT maximum error in calculation of T
c     maximum_iterations,               ! INPUT maximum iterations in calculation of T
c     ncomps,                           ! INPUT total number of components in the state
c     nspecies,                         ! INPUT total number of species in the state
c     state, DIMS(state),               ! INPUT all variables in the state
c     update_for_H, DIMS(update_for_H), ! OUTPUT divergences of the fluxes
c     update_for_Y, DIMS(update_for_Y), ! OUTPUT divergences of the fluxes
c     volume, DIMS(volume),             ! INPUT volumes of the cells
c     xflux_for_H, DIMS(xflux_for_H),   ! OUTPUT x fluxes for enthalpy
c     xflux_for_Y, DIMS(xflux_for_Y),   ! OUTPUT x fluxes for species
c     yflux_for_H, DIMS(yflux_for_H),   ! OUTPUT y fluxes for enthalpy
c     yflux_for_Y, DIMS(yflux_for_Y),   ! OUTPUT y fluxes for species

      subroutine FORT_RK_MULTICOMPONENT
     &     (
     &     domain_lo, domain_hi,
     &     lo, hi,
     &     areax, DIMS(areax),
     &     areay, DIMS(areay),
     &     bc,
     &     dt,
     &     dx,
     &     index_of_firstY,
     &     index_of_lastY,
     &     index_of_rho,
     &     index_of_rhoH,
     &     index_of_T,
     &     maximum_error,
     &     maximum_iterations,
     &     ncomps,
     &     nspecies,
     &     state, DIMS(state),
     &     update_for_H, DIMS(update_for_H),
     &     update_for_Y, DIMS(update_for_Y),
     &     volume, DIMS(volume),
     &     xflux_for_H, DIMS(xflux_for_H),
     &     xflux_for_Y, DIMS(xflux_for_Y),
     &     yflux_for_H, DIMS(yflux_for_H),
     &     yflux_for_Y, DIMS(yflux_for_Y)
     &     )

      implicit none

#include "REAL.H"
#include "BC_TYPES.H"
#include "HEATTRANSFER_F.H"
#include "ChemDriver_F.H"
#include "cdwrk.H"

c     arguments
      integer ncomps
      integer nspecies
      integer domain_lo(BL_SPACEDIM), domain_hi(BL_SPACEDIM)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer DIMDEC(areax)
      integer DIMDEC(areay)
      integer bc(BL_SPACEDIM,2,ncomps)
      integer DIMDEC(state)
      integer DIMDEC(update_for_H)
      integer DIMDEC(update_for_Y)
      integer DIMDEC(xflux_for_H)
      integer DIMDEC(yflux_for_H)
      integer DIMDEC(xflux_for_Y)
      integer DIMDEC(yflux_for_Y)
      integer DIMDEC(volume)
      integer index_of_firstY
      integer index_of_lastY
      integer index_of_rho
      integer index_of_rhoH
      integer index_of_T
      integer maximum_iterations
      REAL_T maximum_error
      REAL_T areax(DIMV(areax))
      REAL_T areay(DIMV(areay))
      REAL_T dt
      REAL_T dx(BL_SPACEDIM)
      REAL_T state(DIMV(state),ncomps)
      REAL_T update_for_H(DIMV(update_for_H))
      REAL_T update_for_Y(DIMV(update_for_Y),nspecies)
      REAL_T xflux_for_H(DIMV(xflux_for_H))
      REAL_T yflux_for_H(DIMV(yflux_for_H))
      REAL_T xflux_for_Y(DIMV(xflux_for_Y),ncomps)
      REAL_T yflux_for_Y(DIMV(yflux_for_Y),ncomps)
      REAL_T volume(DIMV(volume))

c     local variables
      external FORT_TfromHYpt 
      integer FORT_TfromHYpt 
      integer i
      integer iterations
      integer j
      integer n
      integer offset
      logical error
      logical external_check
      logical external_for_T
      logical external_for_Y
      logical internal
      logical internal_check
      REAL_T alpha
      REAL_T ambient_pressure_cgs
      REAL_T ambient_pressure_si
      REAL_T coefficient
      REAL_T contiguous_X
      REAL_T contiguous_Y
      REAL_T cpms
      REAL_T dpdt_factor
      REAL_T dt_over_volume
      REAL_T floored_T
      REAL_T FORT_P1ATMMKS
      REAL_T H
      REAL_T half_over_dx
      REAL_T half_over_dy
      REAL_T hms 
      REAL_T lambda
      REAL_T minus_value
      REAL_T mmw
      REAL_T one_atmosphere_cgs
      REAL_T one_atmosphere_si
      REAL_T one_over_dx
      REAL_T one_over_dy
      REAL_T one_over_rho
      REAL_T one_over_wt
      REAL_T plus_value
      REAL_T residual
      REAL_T rho_cgs
      REAL_T rho_D
      REAL_T rho_Y
      REAL_T RU
      REAL_T RUC
      REAL_T sum
      REAL_T sum_of_fluxes
      REAL_T sum_of_rho_Y
      REAL_T T_on_edge
      REAL_T temperature
      REAL_T Upsilon
      REAL_T X 
      REAL_T Y 
      allocatable contiguous_X
      allocatable contiguous_Y
      allocatable cpms
      allocatable H
      allocatable hms
      allocatable lambda
      allocatable one_over_wt
      allocatable residual
      allocatable rho_D
      allocatable rho_Y
      allocatable temperature
      allocatable Upsilon
      allocatable X
      allocatable Y
      dimension contiguous_X (:)
      dimension contiguous_Y (:)
      dimension cpms (:)
      dimension h (:,:)
      dimension hms (:)
      dimension lambda (:,:)
      dimension one_over_wt (:)
      dimension residual (:)
      dimension rho_D (:,:,:)
      dimension rho_Y (:,:,:)
      dimension temperature (:,:)
      dimension Upsilon (:)
      dimension X (:,:,:)
      dimension Y (:,:,:)

      allocate (contiguous_X (nspecies))
      allocate (contiguous_Y (nspecies))
      allocate (cpms (nspecies))
      allocate (H (DIMV(state)))
      allocate (hms (nspecies))
      allocate (lambda (DIMV(state)))
      allocate (one_over_wt (nspecies))
      allocate (residual (maximum_iterations))
      allocate (rho_D (DIMV(state), nspecies))
      allocate (rho_Y (DIMV(state), nspecies))
      allocate (temperature (DIMV(state)))
      allocate (Upsilon (nspecies))
      allocate (X (DIMV(state), nspecies))
      allocate (Y (DIMV(state), nspecies))

c///
c///  initialization
c///

c     offset of species in the state
      offset = index_of_firstY - 1

c     some constants, at least in Cartesian coordinates
      area_for_x = dx(2)
      area_for_y = dx(1)
      areax_over_dy = dx(1) / dx(2)
      areay_over_dx = dx(2) / dx(1)
      half_over_dx = 0.5d0 / dx(1)
      half_over_dy = 0.5d0 / dx(2)
      one_over_dx = 1.0d0 / dx(1)
      one_over_dy = 1.0d0 / dx(2)

c     get ambient pressure in cgs units
c     one_atmosphere_cgs = 1013250.0 (dyne / cm**2)
c     one_atmosphere_si = 101325.0 (pascal)
      call CKRP (IWRK(ckbi), RWRK(ckbr), universal_gas_constant_cgs, RUC, one_atmosphere_cgs)
      call FORT_GETPAMB (ambient_pressure_si, dpdt_factor)
      one_atmosphere_si = FORT_P1ATMMKS ()
      ambient_pressure_cgs = one_atmosphere_cgs * ambient_pressure_si / one_atmosphere_si

c     copy rho_Y from the state (not strictly necessary but convenient)
c     also, calculate H, T, X, Y from rho, rho H and rho Y in the state
      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            one_over_rho = 1.0d0 / state(i, j, index_of_rho)
            H(i, j) = one_over_rho * state(i, j, index_of_rhoH)
            do n=1,nspecies
               rho_Y(i, j, n) = state(i, j, offset + n)
               Y(i, j, n) = one_over_rho * state(i, j, offset + n)
            end do
            do n=1,nspecies
               contiguous_Y(n) = Y(i, j, n)
            end do
c           use the state value for T as a guess for the Newton iteration
            T(i, j) = state(i, j, index_of_T)
            iterations = FORT_TfromHYpt 
     &           (T(i, j), 
     &           H(i, j), 
     &           contiguous_Y,
     &           maximum_error,
     &           maximum_iterations,
     &           residual)
            CALL CKYTX (contiguous_Y, IWRK(ckbi), RWRK(ckbr), contiguous_X)
            do n=1,nspecies
               X(i, j, n) = contiguous_X(n)
            end do
         end do
      end do

c///
c///  x fluxes
c///

c     set flags for special treatment of the low boundary

      lo_external_for_T = .false.
      lo_external_for_Y = .false.
      lo_internal       = .true.
      if (domain_lo(1) .eq. lo(1)) then
c        reset the flags and check that the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked external then all Y are marked external
c        3) if the Y are marked external then so is T
         lo_external_for_T = EXT_DIR .eq. bc(1, 1, index_of_T)
         lo_external_for_Y = EXT_DIR .eq. bc(1, 1, index_of_firstY)
         lo_external_check = EXT_DIR .eq. bc(1, 1, index_of_firstY)
         lo_internal       = INT_DIR .eq. bc(1, 1, index_of_T)
         lo_internal_check = INT_DIR .eq. bc(1, 1, index_of_T)
         do n=1,nspecies
            lo_external_for_Y = lo_external_for_Y .and. (EXT_DIR .eq. bc(1, 1, offset + n))
            lo_external_check = lo_external_check  .or. (EXT_DIR .eq. bc(1, 1, offset + n))
            lo_internal       = lo_internal       .and. (INT_DIR .eq. bc(1, 1, offset + n))
            lo_internal_check = lo_internal_check  .or. (INT_DIR .eq. bc(1, 1, offset + n))
         end do
         error = .not. (
     &        (lo_external_for_Y .eq. lo_external_check) .and.
     &        (lo_internal .eq. lo_internal_check) .and.
     &        (.not. lo_external_for_Y .or. lo_external_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_DIFFUSION: inconsistent bc array at low x")
      end if

c     set flags for special treatment of the hi boundary

      hi_external_for_T = .false.
      hi_external_for_Y = .false.
      hi_internal       = .true.
      if (domain_hi(1) .eq. hi(1)) then
c        reset the flags and check that the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked external then all Y are marked external
c        3) if the Y are marked external then so is T
         hi_external_for_T = EXT_DIR .eq. bc(1, 2, index_of_T)
         hi_external_for_Y = EXT_DIR .eq. bc(1, 2, index_of_firstY)
         hi_external_check = EXT_DIR .eq. bc(1, 2, index_of_firstY)
         hi_internal       = INT_DIR .eq. bc(1, 2, index_of_T)
         hi_internal_check = INT_DIR .eq. bc(1, 2, index_of_T)
         do n=1,nspecies
            hi_external_for_Y = hi_external_for_Y .and. (EXT_DIR .eq. bc(1, 2, offset + n))
            hi_external_check = hi_external_check  .or. (EXT_DIR .eq. bc(1, 2, offset + n))
            hi_internal       = hi_internal       .and. (INT_DIR .eq. bc(1, 2, offset + n))
            hi_internal_check = hi_internal_check  .or. (INT_DIR .eq. bc(1, 2, offset + n))
         end do
         error = .not. (
     &        (hi_external_for_Y .eq. hi_external_check) .and.
     &        (hi_internal .eq. hi_internal_check) .and.
     &        (.not. hi_external_for_Y .or. hi_external_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_DIFFUSION: inconsistent bc array at high x")
      end if

c     calculate fluxes in the x direction

      do j=lo(2),hi(2)
c        fluxes in x direction are at edges indexed from lo(1) to hi(1)+1
         do i=lo(1),hi(1)+1

c           set flags indicating whether this edge is at the left or right boundary
            at_lo = i .eq. domain_lo(1)
            at_hi = i .eq. domain_hi(1) + 1

c           evaluate the value for T at the edge
c           for external (Dirichlet) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            if (at_lo .and. lo_external_for_T) then
               T_on_edge = state(i, j, index_of_T)
            else if (at_hi .and. hi_external_for_T) then
               T_on_edge = state(i+1, j, index_of_T)
            else
               T_on_edge = 0.5 * (T(i,j) + T(i+1,j))
            end if

c           use a floored T just in case
            T_on_edge = MAX (TMIN_TRANS, T_on_edge)

c           evaluate the values for Y at the edge
c           for external (Dirichlet) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            do n = 1,nspecies
               if (at_lo .and. lo_external_for_Y) then
                  contiguous_Y_on_edge(n) = state(i, j, offset + n)
               else if (at_hi .and. hi_external_for_Y) then
                  contiguous_Y_on_edge(n) = state(i+1, j, offset + n)
               else
                  contiguous_Y_on_edge(n) = 0.5 * (Y(i,j,n) + Y(i+1,j,n)))
               end if
            end do

c           evaluate the values for X at the edge
c           these values are only needed at lo and hi external(Dirichlet) boundaries for Y
            if ((at_lo .and. lo_external_for_Y) .or. (at_hi .and. hi_external_for_Y)) then
               call CKYTX 
     &              (contiguous_Y_on_edge, IWRK(ckbi), RWRK(ckbr), X_on_edge)
            end if

c           evaluate the derivative for T at the edge
            if (at_lo .and. lo_external_for_T) then
               derivative_of_T 
     &              = one_over_dx
     &              * (lo_coefficient_0 * T_on_edge
     &              + lo_coefficient_1 * T(i,j)
     &              + lo_coefficient_2 * T(i+1,j))
            else if (at_hi .and. hi_external_for_T) then
               derivative_of_T 
     &              = one_over_dx
     &              * (hi_coefficient_0 * T(i-2,j)
     &              + hi_coefficient_1 * T(i-1,j)
     &              + hi_coefficient_2 * T_on_edge)
            else
               derivative_of_T = (T(i,j) - T(i-1,j)) * one_over_dx
            end if

c           evaluate the derivatives for X at the edge
            do n = 1,nspecies
               if (at_lo .and. lo_external_for_Y) then
                  dX_over_dx(n)
     &                 = one_over_dx
     &                 * (lo_coefficient_0 * X_on_edge(n)
     &                 + lo_coefficient_1 * X(i,j,n)
     &                 + lo_coefficient_2 * X(i+1,j,n))
               else if (at_hi .and. hi_external_for_Y) then
                  dX_over_dx(n)
     &                 = one_over_dx
     &                 * (hi_coefficient_0 * X(i-2,j,n)
     &                 + hi_coefficient_1 * X(i-1,j,n)
     &                 + hi_coefficient_2 * X_on_edge(n))
               else
                  dX_over_dx(n) = (X(i,j,n) - X(i-1,j,n)) * one_over_dx
               end if
            end do

c           evaluate diffusion coefficients and other things at the edge
            call CKCPMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), cpms_on_edge)
            call EGSPAR 
     &           (T_on_edge, contiguous_X_on_edge, contiguous_Y_on_edge, 
     &           cpms_on_edge, RWRK(egbr), IWRK(egbi))
            call CKMMWY 
     &           (contiguous_Y_on_edge, IWRK(ckbi), RWRK(ckbr), mmw_on_edge)
            call EGSLTDR5 
     &           (T_on_edge, contiguous_Y_on_edge, mmw_on_edge, RWRK(egbr), IWRK(egbi),
     &           lambda_prime_on_edge, rho_theta_on_edge, rho_Y_D_on_edge)
            call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms_on_edge)
c           there does not seem to be a chemkin call to get this ratio from T and Y
            pressure_over_rho_on_edge 
     &           = universal_gas_constant_cgs * T_on_edge / mmw_on_edge

c           assemble the fluxes
            sum_h_xflux_for_Y = 0
            sum_rho_theta_derivative_of_X = 0
            sum_rho_Y_D_derivative_of_X = 0
            do n = 1,nspecies
               sum_rho_Y_D_derivative_of_X = 0
               do ell = 1,nspecies
                  sum_rho_Y_D_derivative_of_X 
     &                 = sum_rho_Y_D_derivative_of_X 
     &                 + rho_Y_D_on_edge(n,ell) 
     &                 * derivative_of_X(ell)
               end do
c              scaling by 0.1 converts expression of mixed units to cgs (check this)
               xflux_for_Y(i,j,n) 
     &              = area_for_x
     &              * (sum_rho_Y_D_derivative_of_X * 1.0d-1
     &              + contiguous_Y_on_edge(n) * rho_theta_on_edge(n) * derivative_of_T / T_on_edge)
               sum_h_xflux_for_Y 
     &              = sum_h_xflux_for_Y 
     &              + hms_on_edge(n) 
     &              * xflux_for_Y(i,j,n)
               sum_rho_theta_derivative_of_X 
     &              = sum_rho_theta_derivative_of_X 
     &              + rho_theta_on_edge(n) 
     &              * derivative_of_X(n)
            end do
c           scaling by 0.1 converts expression of mixed units to cgs (check this)
            xflux_for_H(i,j)
     &           = sum_h_xflux_for_Y
     &           + area_for_x
     &           * (lambda_prime_on_edge * derivative_of_T
     &           + pressure_over_rho_on_edge * sum_rho_theta_derivative_of_X * 1.0d-1)
            end do
         end do

c---- some remarks on scaling

            CALL EGSV1 
     &           (ambient_pressure_cgs, floored_T, contiguous_Y, mmw, RWRK(egbr), 
     &           Upsilon)

            do n=1,nspecies
c              per page 17 in the EGLIB manual, Upsilon has units (cm**2 / s), so
c              rho Upsilon has units (g / cm s).  per equation 12 in the EGLIB manual, 
c              the mixture averaged diffusion coefficient is Upsilon times the 
c              dimensionless mmw / wt(n).  thus rho times the diffusion coefficient 
c              with units (g / (cm s)) = 10 (Kg / (m s)) must be divided by 10 to obtain
c              si units.
               rho_D(i,j,n) = rho_cgs * mmw * one_over_wt(n) * Upsilon(n) * 0.1d0
            end do

            CALL EGSL1 (alpha, floored_T, contiguous_X, RWRK(egbr), plus_value)
            alpha = - 1.0d0
            CALL EGSL1 (alpha, floored_T, contiguous_X, RWRK(egbr), minus_value)
c           per page 15 of the EGLIB manual, these thermal conductivities are returned
c           with units (erg / (cm s K)) = ((J / 10**7) / ((m / 10**2) s K)) =
c           = ((W s) / (10**5 m s K)) = (W / (m K 10**5)).  the si unit of thermal
c           conductivity is (W / (M K)), so these values must be divided by 10**5
c           to convert them to si units.   multiplication by 0.5 averages the values.
            lambda(i,j) = 0.5d-5 * (plus_value + minus_value)

            call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
c           per page 98 of the Chemkin manual, hms has units (erg / g) =
c           = (J / 10**7) / (Kg / 10**3) = (J / Kg) * (1 / 10**4).  so hms
c           must be divided by 10**4 to convert to si units.

c---- end of some remarks on scaling

c     fluxes in x direction are from lo(1) to hi(1)+1
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
c           use a floored temperature just in case
            T_on_edge = MAX (TMIN_TRANS, 0.5 * (T(i,j) + T(i+1,j)))
            do n = 1,nspecies
               contiguous_Y_on_edge(n) = 0.5 * (Y(i,j,n) + Y(i+1,j,n)))
            end do
            call CKCPMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), cpms_on_edge)
            call CKYTX 
     &           (contiguous_Y_on_edge, IWRK(ckbi), RWRK(ckbr), contiguous_X_on_edge)
            call EGSPAR 
     &           (T_on_edge, contiguous_X_on_edge, contiguous_Y_on_edge, 
     &           cpms_on_edge, RWRK(egbr), IWRK(egbi))
            call CKMMWY 
     &           (contiguous_Y_on_edge, IWRK(ckbi), RWRK(ckbr), mmw_on_edge)
            call EGSLTDR5 
     &           (T_on_edge, contiguous_Y_on_edge, mmw_on_edge, RWRK(egbr), IWRK(egbi),
     &           lambda_prime_on_edge, rho_theta_on_edge, rho_Y_D_on_edge)
            call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms_on_edge)
c           there does not seem to be a chemkin call to get this ratio from T and Y
            pressure_over_rho_on_edge 
     &           = universal_gas_constant_cgs * T_on_edge / mmw_on_edge
 
            sum_h_xflux_for_Y = 0
            sum_rho_theta_dX = 0
            sum_rho_Y_D_dX = 0
            do n = 1,nspecies
               sum_rho_Y_D_dX = 0
               do ell = 1,nspecies
                  sum_rho_Y_D_dX 
     &                 = sum_rho_Y_D_dX 
     &                 + rho_Y_D_on_edge(n,ell) 
     &                 * (X(i,j,ell) - X(i-1,j,ell))
               enddo
c              scaling by 0.1 converts expression of mixed units to cgs (check this)
               xflux_for_Y(i,j,n) 
     &              = areay_over_dx 
     &              * (sum_rho_Y_D_dX * 1.0d-1
     &              + contiguous_Y_on_edge(n) * rho_theta_on_edge(n) * (T(i,j) - T(i-1,j)) / T_on_edge)
               sum_h_xflux_for_Y 
     &              = sum_h_xflux_for_Y 
     &              + hms_on_edge(n) 
     &              * xflux_for_Y(i,j,n)
               sum_rho_theta_dX 
     &              = sum_rho_theta_dX 
     &              + rho_theta_on_edge(n) 
     &              * (X(i,j,n) - X(i-1,j,n))
            end do
c           scaling by 0.1 converts expression of mixed units to cgs (check this)
            xflux_for_H(i,j)
     &           = sum_h_xflux_for_Y
     &           + areay_over_dx 
     &           * (lambda_prime_on_edge * (T(i,j) - T(i-1,j))
     &           + pressure_over_rho_on_edge * sum_rho_theta_dX * 1.0d-1)
            end do
         end do

c--------------------------------------------------------------------
c     Note: scalings get cgs quantities from EGLib to MKS (assume He,d come in MKS)
c           sum_rho_Y_D_dX     = zero
c           sum_h_xflux_for_Y    = zero
            sum_rho_theta_dX = zero
            do n = 1,Nspec
c              sum_rho_Y_D_dX = - Ye(i,j,n)*rho_theta_on_edge(n)*d(i,j,Nspec+1)/Te(i,j)
               do L= 1,Nspec
c                 sum_rho_Y_D_dX = sum_rho_Y_D_dX - rhoD(n,L)*d(i,j,L)
               enddo
c              xflux_for_Y(i,j,n) = sum_rho_Y_D_dX*1.d-1
c              sum_h_xflux_for_Y = sum_h_xflux_for_Y + hms_on_edgeH(n)*xflux_for_Y(i,j,n)
c               sum_rho_theta_dX = sum_rho_theta_dX + rho_theta_on_edge(n)*d(i,j,n)
            enddo
            xflux_for_H(i,j) 
c     &           = sum_h_xflux_for_Y
c     &           + (lambda_prime_on_edge*1.d-5)*d(i,j,Nspec+1)
     &           + (RU*Te(i,j)/WW) * (sum_rho_theta_dX * 1.d-1)

c--------------------------------------------------------------------

c           species fluxes
            do n=1,nspecies
               xflux_for_Y(i,j,n) 
     &              = half_over_dx 
     &              * areax(i,j) 
     &              * (rho_D(i-1,j,n) + rho_D(i,j,n)) 
     &              * (Y(i,j,n) - Y(i-1,j,n))
            end do
c           flux correction
            sum_of_fluxes = 0
            sum_of_rho_Y = 0
            do n=1,nspecies
               sum_of_fluxes = sum_of_fluxes + xflux_for_Y(i,j,n)
               sum_of_rho_Y = sum_of_rho_Y + (rho_Y(i-1, j, n) + rho_Y(i, j, n))
            end do
            coefficient = sum_of_fluxes / sum_of_rho_Y
            do n=1,nspecies
               xflux_for_Y(i,j,n) 
     &              = xflux_for_Y(i,j,n) 
     &              - (rho_Y(i-1, j, n) + rho_Y(i, j, n))
     &              * coefficient
            end do

c           enthalpy flux
            T_on_edge = 0.5d0
     &           * (T(i-1,j) 
     &            + T(i,j))
            call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
c           per page 98 of the Chemkin manual, hms has units (erg / g) =
c           = (J / 10**7) / (Kg / 10**3) = (J / Kg) * (1 / 10**4).  so hms
c           must be divided by 10**4 to convert to si units.
            sum = 0
            do n=1,nspecies
               sum = sum + xflux_for_Y(i,j,n) * hms(n) * 1.0d-4
            end do
            xflux_for_H(i,j) 
     &           = sum
     &           + half_over_dx
     &           * areax(i,j) 
     &           * (lambda(i-1,j) + lambda(i,j)) 
     &           * (T(i,j) - T(i-1,j))
         end do
      end do

c     redo the fluxes at low x domain boundary
      if (domain_lo(1) .eq. lo(1)) then
c        check the bc flags
         external_for_T = bc(1, 1, index_of_T) .eq. EXT_DIR
         external_for_Y = bc(1, 1, index_of_firstY) .eq. EXT_DIR
         external_check = bc(1, 1, index_of_firstY) .eq. EXT_DIR
         internal = bc(1, 1, index_of_T) .eq. INT_DIR
         internal_check = bc(1, 1, index_of_T) .eq. INT_DIR
         do n=1,nspecies
            external_for_Y = external_for_Y 
     &           .and. (bc(1, 1, offset + n) .eq. EXT_DIR)
            external_check = external_check 
     &           .or. (bc(1, 1, offset + n) .eq. EXT_DIR)
            internal = internal
     &           .and. (bc(1, 1, offset + n) .eq. INT_DIR)
            internal_check = internal_check 
     &           .or. (bc(1, 1, offset + n) .eq. INT_DIR)
         end do
         error = .not. (
     &        (external_for_Y .eq. external_check) .and.
     &        (internal .eq. internal_check) .and.
     &        (.not. external_for_Y .or. external_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_DIFFUSION: inconsistent bc flags at low x")
         i = lo(1)
         do j=lo(2),hi(2)
            if (external_for_Y) then
c              species fluxes
               do n=1,nspecies
                  xflux_for_Y(i,j,n) 
     &                 = one_over_dx 
     &                 * areax(i,j) 
     &                 * rho_D(i-1,j,n)
     &                 * (- (8.0d0 / 3.0d0) * Y(i-1,j,n)
     &                 + 3 * Y(i,j,n)
     &                 - (1.0d0 / 3.0d0) * Y(i+1,j,n))
               end do
c              flux correction
               sum_of_fluxes = 0
               sum_of_rho_Y = 0
               do n=1,nspecies
                  sum_of_fluxes = sum_of_fluxes + xflux_for_Y(i,j,n)
                  sum_of_rho_Y = sum_of_rho_Y + (rho_Y(i-1, j, n) + rho_Y(i, j, n))
               end do
               coefficient = sum_of_fluxes / sum_of_rho_Y
               do n=1,nspecies
                  xflux_for_Y(i,j,n) 
     &                 = xflux_for_Y(i,j,n) 
     &                 - (rho_Y(i-1, j, n) + rho_Y(i, j, n))
     &                 * coefficient
               end do
            else if (.not. internal) then
               do n=1,nspecies
                  xflux_for_Y(i,j,n) = 0
               end do
            end if
            if (external_for_T) then
c              restore the boundary value
               T(i-1,j) = state(i-1,j,index_of_T)
               T_on_edge = T(i-1,j)
               call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
               sum = 0
               do n=1,nspecies
                  sum = sum + xflux_for_Y(i,j,n) * hms(n) * 1.0d-4
               end do
               xflux_for_H(i,j) 
     &              = sum
     &              + one_over_dx
     &              * areax(i,j) 
     &              * lambda(i-1,j)
     &                 * (- (8.0d0 / 3.0d0) * T(i-1,j)
     &                 + 3 * T(i,j)
     &                 - (1.0d0 / 3.0d0) * T(i+1,j))
            else if (.not. internal) then
               xflux_for_H(i,j) = 0
            end if
         end do

c     redo the fluxes at high x domain boundary
      else if (domain_hi(1) .eq. hi(1)) then
c        check the bc flags
         external_for_T = bc(1, 2, index_of_T) .eq. EXT_DIR
         external_for_Y = bc(1, 2, index_of_firstY) .eq. EXT_DIR
         external_check = bc(1, 2, index_of_firstY) .eq. EXT_DIR
         internal = bc(1, 2, index_of_T) .eq. INT_DIR
         internal_check = bc(1, 2, index_of_T) .eq. INT_DIR
         do n=1,nspecies
            external_for_Y = external_for_Y 
     &           .and. (bc(1, 2, offset + n) .eq. EXT_DIR)
            external_check = external_check 
     &           .or. (bc(1, 2, offset + n) .eq. EXT_DIR)
            internal = internal
     &           .and. (bc(1, 2, offset + n) .eq. INT_DIR)
            internal_check = internal_check 
     &           .or. (bc(1, 2, offset + n) .eq. INT_DIR)
         end do
         error = .not. (
     &        (external_for_Y .eq. external_check) .and.
     &        (internal .eq. internal_check) .and.
     &        (.not. external_for_Y .or. external_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_DIFFUSION: inconsistent bc flags at high x")
         one_over_dx = 1.0d0 / dx(1)
         i = hi(1)
         do j=lo(2),hi(2)
            if (external_for_Y) then
c              species fluxes
               do n=1,nspecies
                  xflux_for_Y(i+1,j,n) 
     &                 = one_over_dx 
     &                 * areax(i+1,j) 
     &                 * rho_D(i+1,j,n)
     &                 * ((1.0d0 / 3.0d0) * Y(i-1,j,n)
     &                 - 3 * Y(i,j,n)
     &                 + (8.0d0 / 3.0d0) * Y(i+1,j,n))
               end do
c              flux correction
               sum_of_fluxes = 0
               sum_of_rho_Y = 0
               do n=1,nspecies
                  sum_of_fluxes = sum_of_fluxes + xflux_for_Y(i+1,j,n)
                  sum_of_rho_Y = sum_of_rho_Y + (rho_Y(i, j, n) + rho_Y(i+1, j, n))
               end do
               coefficient = sum_of_fluxes / sum_of_rho_Y
               do n=1,nspecies
                  xflux_for_Y(i+1,j,n) 
     &                 = xflux_for_Y(i+1,j,n) 
     &                 - (rho_Y(i, j, n) + rho_Y(i+1, j, n))
     &                 * coefficient
               end do
            else if (.not. internal) then
               do n=1,nspecies
                  xflux_for_Y(i+1,j,n) = 0
               end do
            end if
            if (external_for_T) then          
c              restore the boundary value
               T(i+1,j) = state(i+1,j,index_of_T)
               T_on_edge = T(i+1,j)
               call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
               sum = 0
               do n=1,nspecies
                  sum = sum + xflux_for_Y(i+1,j,n) * hms(n) * 1.0d-4
               end do
               xflux_for_H(i+1,j) 
     &              = sum
     &              + one_over_dx
     &              * areax(i+1,j) 
     &              * lambda(i+1,j)
     &              * ((1.0d0 / 3.0d0) * T(i-1,j)
     &              - 3 * T(i,j)
     &              + (8.0d0 / 3.0d0) * T(i+1,j))
            else if (.not. internal) then
               xflux_for_H(i+1,j) = 0
            end if
         end do
      end if

c     fluxes in y direction are from lo(2) to hi(2)+1
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
c           species fluxes
            do n=1,nspecies
               yflux_for_Y(i,j,n) 
     &              = half_over_dy
     &              * areay(i,j) 
     &              * (rho_D(i,j-1,n) + rho_D(i,j,n)) 
     &              * (Y(i,j,n) - Y(i,j-1,n))
            end do
c           flux correction
            sum_of_fluxes = 0
            sum_of_rho_Y = 0
            do n=1,nspecies
               sum_of_fluxes = sum_of_fluxes + yflux_for_Y(i,j,n)
               sum_of_rho_Y = sum_of_rho_Y + (rho_Y(i, j-1, n) + rho_Y(i, j, n))
            end do
            coefficient = sum_of_fluxes / sum_of_rho_Y
            do n=1,nspecies
               yflux_for_Y(i,j,n) 
     &              = yflux_for_Y(i,j,n) 
     &              - (rho_Y(i, j-1, n) + rho_Y(i, j, n))
     &              * coefficient
            end do
c           enthalpy flux
            T_on_edge = 0.5d0
     &           * (T(i,j-1) 
     &            + T(i,j))
            call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
            sum = 0
            do n=1,nspecies
               sum = sum + yflux_for_Y(i,j,n) * hms(n) * 1.0d-4
            end do
            yflux_for_H(i,j) 
     &           = sum
     &           + half_over_dy
     &           * areay(i,j) 
     &           * (lambda(i,j-1) + lambda(i,j)) 
     &           * (T(i,j) - T(i,j-1))
         end do
      end do

c     redo the fluxes at low y domain boundary
      if (domain_lo(2) .eq. lo(2)) then
c        check the bc flags
         external_for_T = bc(2, 1, index_of_T) .eq. EXT_DIR
         external_for_Y = bc(2, 1, index_of_firstY) .eq. EXT_DIR
         external_check = bc(2, 1, index_of_firstY) .eq. EXT_DIR
         internal = bc(2, 1, index_of_T) .eq. INT_DIR
         internal_check = bc(2, 1, index_of_T) .eq. INT_DIR
         do n=1,nspecies
            external_for_Y = external_for_Y 
     &           .and. (bc(2, 1, offset + n) .eq. EXT_DIR)
            external_check = external_check 
     &           .or. (bc(2, 1, offset + n) .eq. EXT_DIR)
            internal = internal
     &           .and. (bc(2, 1, offset + n) .eq. INT_DIR)
            internal_check = internal_check 
     &           .or. (bc(2, 1, offset + n) .eq. INT_DIR)
         end do
         error = .not. (
     &        (external_for_Y .eq. external_check) .and.
     &        (internal .eq. internal_check) .and.
     &        (.not. external_for_Y .or. external_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_DIFFUSION: inconsistent bc flags at low y")
         j = lo(2)
         do i=lo(1),hi(1)
            if (external_for_Y) then
c              species fluxes
               do n=1,nspecies
                  yflux_for_Y(i,j,n) 
     &                 = one_over_dy 
     &                 * areay(i,j) 
     &                 * rho_D(i,j-1,n)
     &                 * (- (8.0d0 / 3.0d0) * Y(i,j-1,n)
     &                 + 3 * Y(i,j,n)
     &                 - (1.0d0 / 3.0d0) * Y(i,j+1,n))
               end do
c              flux correction
               sum_of_fluxes = 0
               sum_of_rho_Y = 0
               do n=1,nspecies
                  sum_of_fluxes = sum_of_fluxes + yflux_for_Y(i,j,n)
                  sum_of_rho_Y = sum_of_rho_Y + (rho_Y(i, j-1, n) + rho_Y(i, j, n))
               end do
               coefficient = sum_of_fluxes / sum_of_rho_Y
               do n=1,nspecies
                  yflux_for_Y(i,j,n) 
     &                 = yflux_for_Y(i,j,n) 
     &                 - (rho_Y(i, j-1, n) + rho_Y(i, j, n))
     &                 * coefficient
               end do
            else if (.not. internal) then
               do n=1,nspecies
                  yflux_for_Y(i,j,n) = 0
               end do
            end if
            if (external_for_T) then          
c              restore the boundary value
               T(i,j-1) = state(i,j-1,index_of_T)
               T_on_edge = T(i,j-1)
               call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
               sum = 0
               do n=1,nspecies
                  sum = sum + yflux_for_Y(i,j,n) * hms(n) * 1.0d-4
               end do
               yflux_for_H(i,j) 
     &              = sum
     &              + one_over_dy
     &              * areay(i,j) 
     &              * lambda(i,j-1)
     &                 * (- (8.0d0 / 3.0d0) * T(i,j-1)
     &                 + 3 * T(i,j)
     &                 - (1.0d0 / 3.0d0) * T(i,j+1))
            else if (.not. internal) then
               yflux_for_H(i,j) = 0
            end if
         end do

c     redo the fluxes at high y domain boundary
      else if (domain_hi(2) .eq. hi(2)) then
c        check the bc flags
         external_for_T = bc(2, 2, index_of_T) .eq. EXT_DIR
         external_for_Y = bc(2, 2, index_of_firstY) .eq. EXT_DIR
         external_check = bc(2, 2, index_of_firstY) .eq. EXT_DIR
         internal = bc(2, 2, index_of_T) .eq. INT_DIR
         internal_check = bc(2, 2, index_of_T) .eq. INT_DIR
         do n=1,nspecies
            external_for_Y = external_for_Y 
     &           .and. (bc(2, 2, offset + n) .eq. EXT_DIR)
            external_check = external_check 
     &           .or. (bc(2, 2, offset + n) .eq. EXT_DIR)
            internal = internal
     &           .and. (bc(2, 2, offset + n) .eq. INT_DIR)
            internal_check = internal_check 
     &           .or. (bc(2, 2, offset + n) .eq. INT_DIR)
         end do
         error = .not. (
     &        (external_for_Y .eq. external_check) .and.
     &        (internal .eq. internal_check) .and.
     &        (.not. external_for_Y .or. external_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_DIFFUSION: inconsistent bc flags at high y")
         one_over_dy = 1.0d0 / dx(2)
         j = hi(2)
         do i=lo(1),hi(1)
            if (external_for_Y) then
c              species fluxes
               do n=1,nspecies
                  yflux_for_Y(i,j+1,n) 
     &                 = one_over_dy 
     &                 * areay(i,j+1) 
     &                 * rho_D(i,j+1,n)
     &                 * ((1.0d0 / 3.0d0) * Y(i,j-1,n)
     &                 - 3 * Y(i,j,n)
     &                 + (8.0d0 / 3.0d0) * Y(i,j+1,n))
               end do
c              flux correction
               sum_of_fluxes = 0
               sum_of_rho_Y = 0
               do n=1,nspecies
                  sum_of_fluxes = sum_of_fluxes + yflux_for_Y(i,j+1,n)
                  sum_of_rho_Y = sum_of_rho_Y + (rho_Y(i, j, n) + rho_Y(i, j+1, n))
               end do
               coefficient = sum_of_fluxes / sum_of_rho_Y
               do n=1,nspecies
                  yflux_for_Y(i,j+1,n) 
     &                 = yflux_for_Y(i,j,n) 
     &                 - (rho_Y(i, j, n) + rho_Y(i, j+1, n))
     &                 * coefficient
               end do
            else if (.not. internal) then
               do n=1,nspecies
                  yflux_for_Y(i,j+1,n) = 0
               end do
            end if
            if (external_for_T) then          
c              restore the boundary value
               T(i,j+1) = state(i,j+1,index_of_T)
               T_on_edge = T(i,j+1)
               call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
               sum = 0
               do n=1,nspecies
                  sum = sum + yflux_for_Y(i,j+1,n) * hms(n) * 1.0d-4
               end do
               yflux_for_H(i,j+1) 
     &              = sum
     &              + one_over_dy
     &              * areay(i,j+1) 
     &              * lambda(i,j+1)
     &              * ((1.0d0 / 3.0d0) * T(i,j-1)
     &              - 3 * T(i,j)
     &              + (8.0d0 / 3.0d0) * T(i,j+1))
            else if (.not. internal) then
               yflux_for_H(i,j+1) = 0
            end if
         end do
      end if

c     volume weighted sum of the fluxes
c     the inclusion of the 1/dx terms leads to failure in the T update in the c++
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            dt_over_volume = dt / volume(i,j)
            update_for_H(i, j) 
     &           = dt_over_volume
     &           * ((xflux_for_H(i+1,j) - xflux_for_H(i,j)) 
     &           + (yflux_for_H(i,j+1) - yflux_for_H(i,j)))
            do n=1,nspecies
               update_for_Y(i, j, n) 
     &              = dt_over_volume
     &              * ((xflux_for_Y(i+1,j,n) - xflux_for_Y(i,j,n))
     &              + (yflux_for_Y(i,j+1,n) - yflux_for_Y(i,j,n)))
            end do
         end do
      end do

      deallocate (contiguous_X)
      deallocate (contiguous_Y)
      deallocate (cpms)
      deallocate (H)
      deallocate (hms)
      deallocate (lambda)
      deallocate (one_over_wt)
      deallocate (residual)
      deallocate (rho_D)
      deallocate (rho_Y)
      deallocate (temperature)
      deallocate (Upsilon)
      deallocate (X)
      deallocate (Y)

      return
      end
