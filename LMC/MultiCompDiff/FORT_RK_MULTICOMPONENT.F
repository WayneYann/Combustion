c#undef BL_LANG_CC
c#ifndef BL_LANG_FORT
c#define BL_LANG_FORT
c#endif

c#define BL_SPACEDIM 2

c///
c///  start of FORT_RK_MULTICOMPONENT
c///

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

#include "ArrayLim.H"
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
      integer ell
      integer i
      integer iterations
      integer j
      integer n
      integer offset
      logical at_hi
      logical at_lo
      logical error
      logical hi_external_check
      logical hi_external_for_T
      logical hi_external_for_Y
      logical hi_internal
      logical hi_internal_check
      logical lo_external_check
      logical lo_external_for_T
      logical lo_external_for_Y
      logical lo_internal
      logical lo_internal_check
      REAL_T alpha
      REAL_T area_for_x
      REAL_T area_for_y
      REAL_T ambient_pressure_cgs
      REAL_T ambient_pressure_si
      REAL_T coefficient
      REAL_T contiguous_X
      REAL_T contiguous_Y
      REAL_T cpms_on_edge
      REAL_T derivative_of_T
      REAL_T derivative_of_X
      REAL_T dpdt_factor
      REAL_T dt_over_volume
      REAL_T floored_T
      REAL_T FORT_P1ATMMKS
      REAL_T H
      REAL_T hi_coefficient_0
      REAL_T hi_coefficient_1
      REAL_T hi_coefficient_2
      REAL_T hms_on_edge 
      REAL_T lambda_prime_on_edge
      REAL_T lo_coefficient_0
      REAL_T lo_coefficient_1
      REAL_T lo_coefficient_2
      REAL_T minus_value
      REAL_T mmw_on_edge
      REAL_T one_atmosphere_cgs
      REAL_T one_atmosphere_si
      REAL_T one_over_dx
      REAL_T one_over_dy
      REAL_T one_over_rho
      REAL_T one_over_wt
      REAL_T plus_value
      REAL_T pressure_over_rho_on_edge
      REAL_T residual
      REAL_T rho_cgs
      REAL_T rho_theta_on_edge
      REAL_T rho_Y_D_on_edge
      REAL_T RU
      REAL_T RUC
      REAL_T sum_h_xflux_for_Y
      REAL_T sum_rho_theta_derivative_of_X
      REAL_T sum_rho_Y_D_derivative_of_X
      REAL_T sum
      REAL_T sum_of_fluxes
      REAL_T sum_of_rho_Y
      REAL_T T
      REAL_T T_on_edge
      REAL_T universal_gas_constant_cgs
      REAL_T Upsilon
      REAL_T X
      REAL_T X_on_edge
      REAL_T Y 
      REAL_T Y_on_edge
      allocatable contiguous_X
      allocatable contiguous_Y
      allocatable cpms_on_edge
      allocatable derivative_of_X
      allocatable H
      allocatable hms_on_edge
      allocatable one_over_wt
      allocatable residual
      allocatable rho_theta_on_edge
      allocatable rho_Y_D_on_edge
      allocatable T
      allocatable Upsilon
      allocatable X
      allocatable X_on_edge
      allocatable Y
      allocatable Y_on_edge
      dimension contiguous_X (:)
      dimension contiguous_Y (:)
      dimension cpms_on_edge (:)
      dimension derivative_of_X (:)
      dimension h (:,:)
      dimension hms_on_edge (:)
      dimension one_over_wt (:)
      dimension residual (:)
      dimension rho_theta_on_edge (:)
      dimension rho_Y_D_on_edge (:,:)
      dimension T (:,:)
      dimension Upsilon (:)
      dimension X (:,:,:)
      dimension X_on_edge (:)
      dimension Y (:,:,:)
      dimension Y_on_edge (:)

      allocate (contiguous_X (nspecies))
      allocate (contiguous_Y (nspecies))
      allocate (cpms_on_edge (nspecies))
      allocate (derivative_of_X (nspecies))
      allocate (H (DIMV(state)))
      allocate (hms_on_edge (nspecies))
      allocate (one_over_wt (nspecies))
      allocate (residual (maximum_iterations))
      allocate (rho_theta_on_edge (nspecies))
      allocate (rho_Y_D_on_edge (nspecies, nspecies))
      allocate (T (DIMV(state)))
      allocate (Upsilon (nspecies))
      allocate (X (DIMV(state), nspecies))
      allocate (X_on_edge (nspecies))
      allocate (Y (DIMV(state), nspecies))
      allocate (Y_on_edge (nspecies))

c///
c///  initialization
c///

c     offset of species in the state
      offset = index_of_firstY - 1

c     some constants, at least in Cartesian coordinates
      area_for_x = dx(2)
      area_for_y = dx(1)
      one_over_dx = 1.0d0 / dx(1)
      one_over_dy = 1.0d0 / dx(2)

c     coefficients for difference formulas at the hi and low edges
      hi_coefficient_0 = + 1.0d0 / 3.0d0
      hi_coefficient_1 = - 3.0d0
      hi_coefficient_2 = + 8.0d0 / 3.0d0
      lo_coefficient_0 = - 8.0d0 / 3.0d0
      lo_coefficient_1 = + 3.0d0
      lo_coefficient_2 = - 1.0d0 / 3.0d0

c     get ambient pressure in cgs units
c     one_atmosphere_cgs = 1013250.0 (dyne / cm**2)
c     one_atmosphere_si = 101325.0 (pascal)
      call CKRP (IWRK(ckbi), RWRK(ckbr), universal_gas_constant_cgs, RUC, one_atmosphere_cgs)
      call FORT_GETPAMB (ambient_pressure_si, dpdt_factor)
      one_atmosphere_si = FORT_P1ATMMKS ()
      ambient_pressure_cgs = one_atmosphere_cgs * ambient_pressure_si / one_atmosphere_si

c     calculate H, T, X, and Y from rho, rho H and rho Y in the state
      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            one_over_rho = 1.0d0 / state(i, j, index_of_rho)
            H(i, j) = one_over_rho * state(i, j, index_of_rhoH)
            do n=1,nspecies
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
                  Y_on_edge(n) = state(i, j, offset + n)
               else if (at_hi .and. hi_external_for_Y) then
                  Y_on_edge(n) = state(i+1, j, offset + n)
               else
                  Y_on_edge(n) = 0.5 * (Y(i,j,n) + Y(i+1,j,n))
               end if
            end do

c           evaluate the values for X at the edge
            call CKYTX 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), X_on_edge)

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
                  derivative_of_X(n)
     &                 = one_over_dx
     &                 * (lo_coefficient_0 * X_on_edge(n)
     &                 + lo_coefficient_1 * X(i,j,n)
     &                 + lo_coefficient_2 * X(i+1,j,n))
               else if (at_hi .and. hi_external_for_Y) then
                  derivative_of_X(n)
     &                 = one_over_dx
     &                 * (hi_coefficient_0 * X(i-2,j,n)
     &                 + hi_coefficient_1 * X(i-1,j,n)
     &                 + hi_coefficient_2 * X_on_edge(n))
               else
                  derivative_of_X(n) = (X(i,j,n) - X(i-1,j,n)) * one_over_dx
               end if
            end do

c           evaluate diffusion coefficients and other things at the edge
            call CKCPMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), cpms_on_edge)
            call EGSPAR 
     &           (T_on_edge, X_on_edge, Y_on_edge, 
     &           cpms_on_edge, RWRK(egbr), IWRK(egbi))
            call CKMMWY 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), mmw_on_edge)
            call EGSLTDR5 
     &           (T_on_edge, Y_on_edge, mmw_on_edge, RWRK(egbr), IWRK(egbi),
     &           lambda_prime_on_edge, rho_theta_on_edge, rho_Y_D_on_edge)
            call CKHMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms_on_edge)
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
     &              + Y_on_edge(n) * rho_theta_on_edge(n) * derivative_of_T / T_on_edge)
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

c            CALL EGSV1 
c     &           (ambient_pressure_cgs, floored_T, contiguous_Y, mmw_on_edge, RWRK(egbr), 
c     &           Upsilon)

c            do n=1,nspecies
c              per page 17 in the EGLIB manual, Upsilon has units (cm**2 / s), so
c              rho Upsilon has units (g / cm s).  per equation 12 in the EGLIB manual, 
c              the mixture averaged diffusion coefficient is Upsilon times the 
c              dimensionless mmw / wt(n).  thus rho times the diffusion coefficient 
c              with units (g / (cm s)) = 10 (Kg / (m s)) must be divided by 10 to obtain
c              si units.
c               rho_D(i,j,n) = rho_cgs * mmw * one_over_wt(n) * Upsilon(n) * 0.1d0
c            end do

c            CALL EGSL1 (alpha, floored_T, contiguous_X, RWRK(egbr), plus_value)
c            alpha = - 1.0d0
c            CALL EGSL1 (alpha, floored_T, contiguous_X, RWRK(egbr), minus_value)
c           per page 15 of the EGLIB manual, these thermal conductivities are returned
c           with units (erg / (cm s K)) = ((J / 10**7) / ((m / 10**2) s K)) =
c           = ((W s) / (10**5 m s K)) = (W / (m K 10**5)).  the si unit of thermal
c           conductivity is (W / (M K)), so these values must be divided by 10**5
c           to convert them to si units.   multiplication by 0.5 averages the values.
c            lambda(i,j) = 0.5d-5 * (plus_value + minus_value)

c            call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
c           per page 98 of the Chemkin manual, hms has units (erg / g) =
c           = (J / 10**7) / (Kg / 10**3) = (J / Kg) * (1 / 10**4).  so hms
c           must be divided by 10**4 to convert to si units.

c---- end of some remarks on scaling

c///
c///  form the updates weighted by dt / vol
c///

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

c///
c///  release temporary arrays (not strictly necessary)
c///

      deallocate (contiguous_X)
      deallocate (contiguous_Y)
      deallocate (cpms_on_edge)
      deallocate (derivative_of_X)
      deallocate (H)
      deallocate (hms_on_edge)
      deallocate (one_over_wt)
      deallocate (residual)
      deallocate (rho_theta_on_edge)
      deallocate (T)
      deallocate (Upsilon)
      deallocate (X)
      deallocate (X_on_edge)
      deallocate (Y)
      deallocate (Y_on_edge)

c///
c///  end of FORT_RK_MULTICOMPONENT
c///

      return
      end
