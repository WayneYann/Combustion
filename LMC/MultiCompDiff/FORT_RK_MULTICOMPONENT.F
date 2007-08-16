c#undef BL_LANG_CC
c#ifndef BL_LANG_FORT
c#define BL_LANG_FORT
c#endif

c#define BL_SPACEDIM 2

c///
c///  start of FORT_RK_MULTICOMPONENT
c///

c     evaluate the "operator" (extensive fluxes and their divergences) for the 
c     Runge-Kutta implementation of the diffusion update using the multicomponent
c     formulation with Soret and Dufour effects

c     the divergences are scaled by 1/vol so they do approximate the analytic
c     divergence.  when used in the Runge-Kutta formula their additional scaling,
c     given by scale_factor, should be set to - dt.  the minus occurs because
c     the standard ODE is written y_prime = f(y), that is, the divergences are
c     moved to the opposite side of the equation from the time derivative.  when
c     used in getViscTerms the scale_factor should be set to -1.

c     to do:
c     1) must adjust diffusion coefficients by thickness factors
c     2) unsure of LMC support for r-z coordinates

c     arguments are alphabetical, mostly:
c
c     domain_lo, domain_hi,             ! INPUT limits of valid region of the domain
c     lo, hi,                           ! INPUT limits of valid region of the box
c   * areax, DIMS(areax),               ! INPUT areas of the faces perendicular to x axis
c   * areay, DIMS(areay),               ! INPUT areas of the faces perpendicular to y axis
c     bc,                               ! INPUT boundary condition array for all comps
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
c     scale_factor,                     ! INPUT scale_factor
c     state, DIMS(state),               ! INPUT all variables in the state
c   * volume, DIMS(volume),             ! INPUT volumes of the cells
c     div_of_flux_for_H, DIMS(div_of_flux_for_H), ! OUTPUT divergence of the flux for rho H
c     div_of_flux_for_Y, DIMS(div_of_flux_for_Y), ! OUTPUT divergences of the fluxes for rho Y
c     xflux_for_H, DIMS(xflux_for_H),   ! OUTPUT extensive x fluxes for rho H
c     xflux_for_Y, DIMS(xflux_for_Y),   ! OUTPUT extensive x fluxes for rho Y
c     yflux_for_H, DIMS(yflux_for_H),   ! OUTPUT extensive y fluxes for rho H
c     yflux_for_Y, DIMS(yflux_for_Y),   ! OUTPUT extensive y fluxes for rho Y
c
c     * these arguments are not used

      subroutine FORT_RK_MULTICOMPONENT
     &     (
     &     domain_lo, domain_hi,
     &     lo, hi,
     &     areax, DIMS(areax),
     &     areay, DIMS(areay),
     &     bc,
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
     &     scale_factor,
     &     state, DIMS(state),
     &     volume, DIMS(volume),
     &     div_of_flux_for_H, DIMS(div_of_flux_for_H),
     &     div_of_flux_for_Y, DIMS(div_of_flux_for_Y),
     &     xflux_for_H, DIMS(xflux_for_H),
     &     xflux_for_Y, DIMS(xflux_for_Y),
     &     yflux_for_H, DIMS(yflux_for_H),
     &     yflux_for_Y, DIMS(yflux_for_Y)
     &     )

      implicit none

#include "BC_TYPES.H"
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
      integer DIMDEC(div_of_flux_for_H)
      integer DIMDEC(div_of_flux_for_Y)
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
      REAL_T areax(DIMV(areax))
      REAL_T areay(DIMV(areay))
      REAL_T div_of_flux_for_H(DIMV(div_of_flux_for_H))
      REAL_T div_of_flux_for_Y(DIMV(div_of_flux_for_Y),nspecies)
      REAL_T dx(BL_SPACEDIM)
      REAL_T maximum_error
      REAL_T scale_factor
      REAL_T state(DIMV(state),ncomps)
      REAL_T volume(DIMV(volume))
      REAL_T xflux_for_H(DIMV(xflux_for_H))
      REAL_T yflux_for_H(DIMV(yflux_for_H))
      REAL_T xflux_for_Y(DIMV(xflux_for_Y),ncomps)
      REAL_T yflux_for_Y(DIMV(yflux_for_Y),ncomps)

c     local variables
      external FORT_TfromHYpt 
      integer FORT_TfromHYpt
      integer ell
      integer i
      integer iterations
      integer j
      integer n
      integer ni
      integer nj
      integer offset
      logical at_hi
      logical at_lo
      logical error
      logical hi_dirichlet_check
      logical hi_dirichlet_for_T
      logical hi_dirichlet_for_Y
      logical hi_internal
      logical hi_internal_check
      logical lo_dirichlet_check
      logical lo_dirichlet_for_T
      logical lo_dirichlet_for_Y
      logical lo_internal
      logical lo_internal_check
      REAL_T area_for_x
      REAL_T area_for_y
      REAL_T ambient_pressure_cgs
      REAL_T ambient_pressure_si
      REAL_T contiguous_X
      REAL_T contiguous_Y
      REAL_T cpms_on_edge
      REAL_T derivative_of_T
      REAL_T derivative_of_X
      REAL_T dpdt_factor
      REAL_T floored_T
      REAL_T FORT_P1ATMMKS
      REAL_T H
      REAL_T hi_coefficient_0
      REAL_T hi_coefficient_1
      REAL_T hi_coefficient_2
      REAL_T hms_on_edge_cgs
      REAL_T hms_on_edge_si
      REAL_T lambda_prime_on_edge_cgs
      REAL_T lambda_prime_on_edge_si
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
      REAL_T plus_value
      REAL_T pressure_over_rho_on_edge_cgs
      REAL_T pressure_over_rho_on_edge_si
      REAL_T residual
      REAL_T rho_theta_on_edge_cgs
      REAL_T rho_theta_on_edge_si
      REAL_T rho_Y_D_on_edge_cgs
      REAL_T rho_Y_D_on_edge_si
      REAL_T RUC
      REAL_T scale_factor_over_volume
      REAL_T sum_h_xflux_for_Y
      REAL_T sum_h_yflux_for_Y
      REAL_T sum_rho_theta_derivative_of_X
      REAL_T sum_rho_Y_D_derivative_of_X
      REAL_T T
      REAL_T T_on_edge
      REAL_T universal_gas_constant_cgs
      REAL_T X
      REAL_T X_on_edge
      REAL_T Y 
      REAL_T Y_on_edge
      allocatable contiguous_X
      allocatable contiguous_Y
      allocatable cpms_on_edge
      allocatable derivative_of_X
      allocatable hms_on_edge_cgs
      allocatable hms_on_edge_si
      allocatable residual
      allocatable rho_theta_on_edge_cgs
      allocatable rho_theta_on_edge_si
      allocatable rho_Y_D_on_edge_cgs
      allocatable rho_Y_D_on_edge_si
      allocatable T
      allocatable X
      allocatable X_on_edge
      allocatable Y
      allocatable Y_on_edge
      dimension contiguous_X (:)
      dimension contiguous_Y (:)
      dimension cpms_on_edge (:)
      dimension derivative_of_X (:)
      dimension hms_on_edge_cgs (:)
      dimension hms_on_edge_si (:)
      dimension residual (:)
      dimension rho_theta_on_edge_cgs (:)
      dimension rho_theta_on_edge_si (:)
      dimension rho_Y_D_on_edge_cgs (:,:)
      dimension rho_Y_D_on_edge_si (:,:)
      dimension T (:,:)
      dimension X (:,:,:)
      dimension X_on_edge (:)
      dimension Y (:,:,:)
      dimension Y_on_edge (:)

c     debug variables
      integer idx
      integer jdx
      logical debug
      logical debug2
c     debug prints patches
      debug = .false.
      debug2 = .false.
      idx = hi(1)
      jdx = hi(2)
      idx = 3
      jdx = 127

      allocate (contiguous_X (nspecies))
      allocate (contiguous_Y (nspecies))
      allocate (cpms_on_edge (nspecies))
      allocate (derivative_of_X (nspecies))
      allocate (hms_on_edge_cgs (nspecies))
      allocate (hms_on_edge_si (nspecies))
      allocate (residual (maximum_iterations))
      allocate (rho_theta_on_edge_cgs (nspecies))
      allocate (rho_theta_on_edge_si (nspecies))
      allocate (rho_Y_D_on_edge_cgs (nspecies, nspecies))
      allocate (rho_Y_D_on_edge_si (nspecies, nspecies))
      allocate (T (DIMV(state)))
      allocate (X (DIMV(state), nspecies))
      allocate (X_on_edge (nspecies))
      allocate (Y (DIMV(state), nspecies))
      allocate (Y_on_edge (nspecies))

c///
c///  initialization
c///

      if (debug) print *, "JFG: at top of FORT_RK_MULTICOMPONENT"

c     offset of species in the state
      offset = index_of_firstY - 1

c     some constants, at least in Cartesian coordinates
      area_for_x = dx(2)
      area_for_y = dx(1)
      one_over_dx = 1.0d0 / dx(1)
      one_over_dy = 1.0d0 / dx(2)
      scale_factor_over_volume = scale_factor / (dx(1) * dx(2))

c     coefficients for difference formulas at the hi and low edges
c     note still must divide by dx or dy to get the approximate derivative
      hi_coefficient_0 = + 1.0d0 / 3.0d0
      hi_coefficient_1 = - 3.0d0
      hi_coefficient_2 = + 8.0d0 / 3.0d0
      lo_coefficient_0 = - 8.0d0 / 3.0d0
      lo_coefficient_1 = + 3.0d0
      lo_coefficient_2 = - 1.0d0 / 3.0d0

c     get ambient pressure
c     one_atmosphere_cgs = 1013250.0 (dyne / cm**2) = 10.0 (Pascal)
c     one_atmosphere_si = 101325.0 (Pascal)
      call CKRP (IWRK(ckbi), RWRK(ckbr), universal_gas_constant_cgs, RUC, one_atmosphere_cgs)
      call FORT_GETPAMB (ambient_pressure_si, dpdt_factor)
      one_atmosphere_si = FORT_P1ATMMKS ()
      ambient_pressure_cgs = one_atmosphere_cgs * ambient_pressure_si / one_atmosphere_si

c     calculate T, X, and Y from rho, rho H and rho Y in the state
      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            one_over_rho = 1.0d0 / state(i, j, index_of_rho)
            H = one_over_rho * state(i, j, index_of_rhoH)
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
     &           H, 
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

      if (debug) call print_patch_n ("T", T, DIMS(state), 1, idx, jdx)
      if (debug) call print_patch_n ("X", X, DIMS(state), nspecies, idx, jdx)

c///
c///  x fluxes
c///

      if (debug2) then
         print *
         print *, "////////// x fluxes //////////"
         print *
      end if

c     set flags for special treatment of the lo boundary

      lo_dirichlet_for_T = .false.
      lo_dirichlet_for_Y = .false.
      lo_internal        = .true.
      if (domain_lo(1) .eq. lo(1)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         lo_dirichlet_for_T = EXT_DIR .eq. bc(1, 1, index_of_T)
         lo_dirichlet_for_Y = EXT_DIR .eq. bc(1, 1, index_of_firstY)
         lo_dirichlet_check = EXT_DIR .eq. bc(1, 1, index_of_firstY)
         lo_internal        = INT_DIR .eq. bc(1, 1, index_of_T)
         lo_internal_check  = INT_DIR .eq. bc(1, 1, index_of_T)
         do n=1,nspecies
            lo_dirichlet_for_Y = lo_dirichlet_for_Y .and. (EXT_DIR .eq. bc(1, 1, offset + n))
            lo_dirichlet_check = lo_dirichlet_check  .or. (EXT_DIR .eq. bc(1, 1, offset + n))
            lo_internal        = lo_internal        .and. (INT_DIR .eq. bc(1, 1, offset + n))
            lo_internal_check  = lo_internal_check   .or. (INT_DIR .eq. bc(1, 1, offset + n))
         end do
         error = .not. (
     &        (lo_dirichlet_for_Y .eqv. lo_dirichlet_check) .and.
     &        (lo_internal .eqv. lo_internal_check) .and.
     &        (.not. lo_dirichlet_for_Y .or. lo_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_DIFFUSION: inconsistent bc array at low x")
      end if

c     set flags for special treatment of the hi boundary

      hi_dirichlet_for_T = .false.
      hi_dirichlet_for_Y = .false.
      hi_internal        = .true.
      if (domain_hi(1) .eq. hi(1)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         hi_dirichlet_for_T = EXT_DIR .eq. bc(1, 2, index_of_T)
         hi_dirichlet_for_Y = EXT_DIR .eq. bc(1, 2, index_of_firstY)
         hi_dirichlet_check = EXT_DIR .eq. bc(1, 2, index_of_firstY)
         hi_internal        = INT_DIR .eq. bc(1, 2, index_of_T)
         hi_internal_check  = INT_DIR .eq. bc(1, 2, index_of_T)
         do n=1,nspecies
            hi_dirichlet_for_Y = hi_dirichlet_for_Y .and. (EXT_DIR .eq. bc(1, 2, offset + n))
            hi_dirichlet_check = hi_dirichlet_check  .or. (EXT_DIR .eq. bc(1, 2, offset + n))
            hi_internal        = hi_internal        .and. (INT_DIR .eq. bc(1, 2, offset + n))
            hi_internal_check  = hi_internal_check   .or. (INT_DIR .eq. bc(1, 2, offset + n))
         end do
         error = .not. (
     &        (hi_dirichlet_for_Y .eqv. hi_dirichlet_check) .and.
     &        (hi_internal .eqv. hi_internal_check) .and.
     &        (.not. hi_dirichlet_for_Y .or. hi_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_DIFFUSION: inconsistent bc array at high x")
      end if

c     calculate fluxes in the x direction

      do j=lo(2),hi(2)
c        fluxes in x direction are at edges indexed from lo(1) to hi(1)+1
         do i=lo(1),hi(1)+1

c           set flags indicating whether the edge is at the left or right boundary
            at_lo = i .eq. domain_lo(1)
            at_hi = i .eq. domain_hi(1) + 1

c           evaluate the value for T on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            if (at_lo .and. lo_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i-1, j, index_of_T)
            else if (at_hi .and. hi_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i, j, index_of_T)
            else
               T_on_edge = 0.5d0 * (T(i-1,j) + T(i,j))
            end if

c           use a floored T just in case
            T_on_edge = MAX (TMIN_TRANS, T_on_edge)

c           evaluate the values for Y on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            do n = 1,nspecies
               if (at_lo .and. lo_dirichlet_for_Y) then
                  Y_on_edge(n) = Y(i-1, j, n)
               else if (at_hi .and. hi_dirichlet_for_Y) then
                  Y_on_edge(n) = Y(i, j, n)
               else
                  Y_on_edge(n) = 0.5d0 * (Y(i-1,j,n) + Y(i,j,n))
               end if
            end do

c           evaluate the values for X on the edge
            call CKYTX 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), X_on_edge)

            if (debug2 .and. i .eq. idx .and. j .eq. jdx) then
               print *, "T_on_edge", T_on_edge
               print *, "Y_on_edge", Y_on_edge
               print *, "X_on_edge", X_on_edge
            end if

c           evaluate the derivative for T on the edge
            if (at_lo .and. .not. lo_internal) then
               if (lo_dirichlet_for_T) then
                  derivative_of_T 
     &                 = one_over_dx
     &                 * (lo_coefficient_0 * T_on_edge
     &                 + lo_coefficient_1 * T(i,j)
     &                 + lo_coefficient_2 * T(i+1,j))
               else
                  derivative_of_T = 0.0d0
               end if
            else if (at_hi .and. .not. hi_internal) then
               if (hi_dirichlet_for_T) then
                  derivative_of_T 
     &                 = one_over_dx
     &                 * (hi_coefficient_0 * T(i-2,j)
     &                 + hi_coefficient_1 * T(i-1,j)
     &                 + hi_coefficient_2 * T_on_edge)
               else
                  derivative_of_T = 0.0d0
               end if
            else
               derivative_of_T = (T(i,j) - T(i-1,j)) * one_over_dx
            end if

c           evaluate the derivatives for X on the edge
            do n = 1,nspecies
               if (at_lo .and. .not. lo_internal) then
                  if (lo_dirichlet_for_Y) then
                     derivative_of_X(n)
     &                    = one_over_dx
     &                    * (lo_coefficient_0 * X_on_edge(n)
     &                    + lo_coefficient_1 * X(i,j,n)
     &                    + lo_coefficient_2 * X(i+1,j,n))
                  else
                     derivative_of_X(n) = 0.0d0
                  end if
               else if (at_hi .and. .not. hi_internal) then
                  if (hi_dirichlet_for_Y) then
                     derivative_of_X(n)
     &                    = one_over_dx
     &                    * (hi_coefficient_0 * X(i-2,j,n)
     &                    + hi_coefficient_1 * X(i-1,j,n)
     &                    + hi_coefficient_2 * X_on_edge(n))
                  else
                     derivative_of_X(n) = 0.0d0
                  end if
               else
                  derivative_of_X(n) = (X(i,j,n) - X(i-1,j,n)) * one_over_dx
               end if
            end do

            if (debug2 .and. i .eq. idx .and. j .eq. jdx) then
               print *, "derivative_of_T", derivative_of_T
               print *, "derivative_of_X", derivative_of_X
            end if

c           evaluate diffusion coefficients and other things on the edge
            call CKCPMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), cpms_on_edge)
            call EGSPAR 
     &           (T_on_edge, X_on_edge, Y_on_edge, 
     &           cpms_on_edge, RWRK(egbr), IWRK(egbi))
            call CKMMWY 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), mmw_on_edge)
c           SUBROUTINE EGSLTDR5 (T, Y, WW, WEG, IWEG, PTC, THETA, D)
            call EGSLTDR5
     &           (T_on_edge, Y_on_edge, mmw_on_edge, RWRK(egbr), IWRK(egbi),
     &           lambda_prime_on_edge_cgs, rho_theta_on_edge_cgs, rho_Y_D_on_edge_cgs)
            call CKHMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms_on_edge_cgs)
c           there does not seem to be a chemkin call to get this ratio from T and Y
            pressure_over_rho_on_edge_cgs 
     &           = universal_gas_constant_cgs * T_on_edge / mmw_on_edge

            if (debug2 .and. i .eq. idx .and. j .eq. jdx) then
               print *, "cpms_on_edge", cpms_on_edge
               print *, "mmw_on_edge", mmw_on_edge
               print *, "hms_on_edge_cgs", hms_on_edge_cgs
               print *, "lambda_prime_on_edge_cgs", lambda_prime_on_edge_cgs
               print *, "rho_theta_on_edge_cgs", rho_theta_on_edge_cgs
               print *, "rho_Y_D_on_edge_cgs", rho_Y_D_on_edge_cgs
            end if

c           convert the library values from cgs units to si units

c           per page 98 of the chemkin manual, the enthalpies for the species,
c           hms_on_edge_cgs, have units
c
c                            erg / g
c                 =          (J / 10**7) / (Kg / 10**3)
c                 = 1.0d-4 * J / Kg
c
c           the si unit of thermal conductivity is J / Kg, so the EGLIB values must
c           be multiplied by 1.0d-4 to convert them to si units.
            do n=1,nspecies
               hms_on_edge_si(n) = 1.0d-4 * hms_on_edge_cgs(n)
            end do

c           per page 23 of the EGLIB manual, the partial thermal conductivity, 
c           lambda_prime_on_edge_cgs, has units
c
c                            erg / (cm s K)
c                 =          (J / 10**7) / ((m / 10**2) s K)
c                 = 1.0d-5 * ((W s) / (m s K))
c                 = 1.0d-5 * W / (m K).
c
c           the si unit of thermal conductivity is W / (m K), so the EGLIB value must
c           be multiplied by 1.0d-5 to convert it to si units.
            lambda_prime_on_edge_si = 1.0d-5 * lambda_prime_on_edge_cgs

c           per page 23 in the EGLIB manual, the rescaled thermal diffusion vector,
c           rho_theta_on_edge_cgs, has units
c
c                            g / (cm s)
c                 =          (Kg / 10**3) / ((m / 10**2) s)
c                 = 1.0d-1 * Kg / (m s)
c
c           so the EGLIB values must be multiplied by 1.0d-1 to convert them to si units.
            do n=1,nspecies
               rho_theta_on_edge_si(n) = 1.0d-1 * rho_theta_on_edge_cgs(n)
            end do

c           per page 23 in the EGLIB manual, the rescaled flux diffusion matrix,
c           rho_Y_D_on_edge_cgs, has units
c
c                            g / (cm s)
c                 =          (Kg / 10**3) / ((m / 10**2) s)
c                 = 1.0d-1 * Kg / (m s)
c
c           so the EGLIB values must be multiplied by 1.0d-1 to convert them to si units.
            do ni=1,nspecies
               do nj=1,nspecies
                  rho_Y_D_on_edge_si(ni,nj) = 1.0d-1 * rho_Y_D_on_edge_cgs(ni,nj)
               end do
            end do

c           the cgs unit for pressure is the barye which equals 0.1 Pascal.
c           the si unit for pressure is the Pascal which equals 1 Kg m s**2.

c           the pressure and density ratio, pressure_over_rho_on_edge_cgs, has units
c
c                            barye / (g / cm**3)
c                 =          (1.0d-1 Pa) / ((1.0d-3 Kg) / (1.0d-2 m)**3)
c                 = 1.0d-4 * Pa / (Kg / m**2)
c
c           so the cgs value must be multiplied by 1.0d-4 to convert it to si units.
            pressure_over_rho_on_edge_si = 1.0d-4 * pressure_over_rho_on_edge_cgs 

            if (debug2 .and. i .eq. idx .and. j .eq. jdx) then
               print *, "hms_on_edge_si", hms_on_edge_si
               print *, "lambda_prime_on_edge_si", lambda_prime_on_edge_si
               print *, "rho_theta_on_edge_si", rho_theta_on_edge_si
               print *, "rho_Y_D_on_edge_si", rho_Y_D_on_edge_si
            end if

c           assemble the fluxes
            sum_h_xflux_for_Y = 0
            sum_rho_theta_derivative_of_X = 0
            do n = 1,nspecies
               sum_rho_Y_D_derivative_of_X = 0
               do ell = 1,nspecies
                  sum_rho_Y_D_derivative_of_X 
     &                 = sum_rho_Y_D_derivative_of_X 
     &                 + rho_Y_D_on_edge_si(n,ell) 
     &                 * derivative_of_X(ell)
               end do
               xflux_for_Y(i,j,n) 
     &              = - area_for_x
     &              * (sum_rho_Y_D_derivative_of_X
     &              + Y_on_edge(n) * rho_theta_on_edge_si(n) * derivative_of_T / T_on_edge)
               sum_h_xflux_for_Y 
     &              = sum_h_xflux_for_Y 
     &              + hms_on_edge_si(n) 
     &              * xflux_for_Y(i,j,n)
               sum_rho_theta_derivative_of_X 
     &              = sum_rho_theta_derivative_of_X 
     &              + rho_theta_on_edge_si(n)
     &              * derivative_of_X(n)
            end do
            xflux_for_H(i,j)
     &           = sum_h_xflux_for_Y
     &           - area_for_x
     &           * (lambda_prime_on_edge_si * derivative_of_T
     &           + pressure_over_rho_on_edge_si * sum_rho_theta_derivative_of_X)
            end do
         end do

      if (debug) call print_patch_n ("xflux_for_H", xflux_for_H, DIMS(xflux_for_H), 1, idx, jdx)
      if (debug) call print_patch_n ("xflux_for_Y", xflux_for_Y, DIMS(xflux_for_Y), nspecies, idx, jdx)

c///
c///  y fluxes
c///

      if (debug2) then
         print *
         print *, "////////// y fluxes //////////"
         print *
      end if

c     set flags for special treatment of the lo boundary

      lo_dirichlet_for_T = .false.
      lo_dirichlet_for_Y = .false.
      lo_internal        = .true.
      if (domain_lo(2) .eq. lo(2)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         lo_dirichlet_for_T = EXT_DIR .eq. bc(2, 1, index_of_T)
         lo_dirichlet_for_Y = EXT_DIR .eq. bc(2, 1, index_of_firstY)
         lo_dirichlet_check = EXT_DIR .eq. bc(2, 1, index_of_firstY)
         lo_internal        = INT_DIR .eq. bc(2, 1, index_of_T)
         lo_internal_check  = INT_DIR .eq. bc(2, 1, index_of_T)
         do n=1,nspecies
            lo_dirichlet_for_Y = lo_dirichlet_for_Y .and. (EXT_DIR .eq. bc(2, 1, offset + n))
            lo_dirichlet_check = lo_dirichlet_check  .or. (EXT_DIR .eq. bc(2, 1, offset + n))
            lo_internal        = lo_internal        .and. (INT_DIR .eq. bc(2, 1, offset + n))
            lo_internal_check  = lo_internal_check   .or. (INT_DIR .eq. bc(2, 1, offset + n))
         end do
         error = .not. (
     &        (lo_dirichlet_for_Y .eqv. lo_dirichlet_check) .and.
     &        (lo_internal .eqv. lo_internal_check) .and.
     &        (.not. lo_dirichlet_for_Y .or. lo_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_DIFFUSION: inconsistent bc array at low y")
      end if

c     set flags for special treatment of the hi boundary

      hi_dirichlet_for_T = .false.
      hi_dirichlet_for_Y = .false.
      hi_internal        = .true.
      if (domain_hi(2) .eq. hi(2)) then
c        set flags and check whether the bc array is consistent for T and Y.
c        consistent means:
c        1) if T or one Y is marked internal then T and all Y are marked internal
c        2) if one Y is marked dirichlet then all Y are marked dirichlet
c        3) if the Y are marked dirichlet then so is T
         hi_dirichlet_for_T = EXT_DIR .eq. bc(2, 2, index_of_T)
         hi_dirichlet_for_Y = EXT_DIR .eq. bc(2, 2, index_of_firstY)
         hi_dirichlet_check = EXT_DIR .eq. bc(2, 2, index_of_firstY)
         hi_internal        = INT_DIR .eq. bc(2, 2, index_of_T)
         hi_internal_check  = INT_DIR .eq. bc(2, 2, index_of_T)
         do n=1,nspecies
            hi_dirichlet_for_Y = hi_dirichlet_for_Y .and. (EXT_DIR .eq. bc(2, 2, offset + n))
            hi_dirichlet_check = hi_dirichlet_check  .or. (EXT_DIR .eq. bc(2, 2, offset + n))
            hi_internal        = hi_internal        .and. (INT_DIR .eq. bc(2, 2, offset + n))
            hi_internal_check  = hi_internal_check   .or. (INT_DIR .eq. bc(2, 2, offset + n))
         end do
         error = .not. (
     &        (hi_dirichlet_for_Y .eqv. hi_dirichlet_check) .and.
     &        (hi_internal .eqv. hi_internal_check) .and.
     &        (.not. hi_dirichlet_for_Y .or. hi_dirichlet_for_T))
         if (error) call bl_abort 
     &        ("FORT_RK_DIFFUSION: inconsistent bc array at high y")
      end if

c     calculate fluxes in the y direction

c     fluxes in y direction are at edges indexed from lo(2) to hi(2)+1
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)

c           set flags indicating whether the edge is at the left or right boundary
            at_lo = j .eq. domain_lo(2)
            at_hi = j .eq. domain_hi(2) + 1

c           evaluate the value for T on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            if (at_lo .and. lo_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i, j-1, index_of_T)
            else if (at_hi .and. hi_dirichlet_for_T) then
c              take the boundary value in state rather than the computed value in T
               T_on_edge = state(i, j, index_of_T)
            else
               T_on_edge = 0.5d0 * (T(i,j-1) + T(i,j))
            end if

c           use a floored T just in case
            T_on_edge = MAX (TMIN_TRANS, T_on_edge)

c           evaluate the values for Y on the edge
c           for Dirichlet (external) boundaries the fill patch iterator has
c           placed the boundary value in the ghost cell adjacent to the edge
            do n = 1,nspecies
               if (at_lo .and. lo_dirichlet_for_Y) then
                  Y_on_edge(n) = Y(i, j-1, n)
               else if (at_hi .and. hi_dirichlet_for_Y) then
                  Y_on_edge(n) = Y(i, j, n)
               else
                  Y_on_edge(n) = 0.5d0 * (Y(i,j-1,n) + Y(i,j,n))
               end if
            end do

c           evaluate the values for X on the edge
            call CKYTX
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), X_on_edge)

            if (debug2 .and. i .eq. idx .and. j .eq. jdx) then
               print *, "T_on_edge", T_on_edge
               print *, "Y_on_edge", Y_on_edge
               print *, "X_on_edge", X_on_edge
            end if

c           evaluate the derivative for T on the edge
            if (at_lo .and. .not. lo_internal) then
               if (lo_dirichlet_for_T) then
                  derivative_of_T
     &                 = one_over_dy
     &                 * (lo_coefficient_0 * T_on_edge
     &                 + lo_coefficient_1 * T(i,j)
     &                 + lo_coefficient_2 * T(i,j+1))
               else
                  derivative_of_T = 0
               end if
            else if (at_hi .and. .not. hi_internal) then
               if (hi_dirichlet_for_T) then
                  derivative_of_T
     &                 = one_over_dy
     &                 * (hi_coefficient_0 * T(i,j-2)
     &                 + hi_coefficient_1 * T(i,j-1)
     &                 + hi_coefficient_2 * T_on_edge)
               else
                  derivative_of_T = 0
               end if
            else
               derivative_of_T = (T(i,j) - T(i,j-1)) * one_over_dy
            end if

c           evaluate the derivatives for X on the edge
            do n = 1,nspecies
               if (at_lo .and. .not. lo_internal) then
                  if (lo_dirichlet_for_Y) then
                     derivative_of_X(n)
     &                    = one_over_dy
     &                    * (lo_coefficient_0 * X_on_edge(n)
     &                    + lo_coefficient_1 * X(i,j,n)
     &                    + lo_coefficient_2 * X(i,j+1,n))
                  else
                     derivative_of_X(n) = 0.0d0
                  end if                     
               else if (at_hi .and. .not. hi_internal) then
                  if (hi_dirichlet_for_Y) then
                     derivative_of_X(n)
     &                    = one_over_dy
     &                    * (hi_coefficient_0 * X(i,j-2,n)
     &                    + hi_coefficient_1 * X(i,j-1,n)
     &                    + hi_coefficient_2 * X_on_edge(n))
                  else
                     derivative_of_X(n) = 0.0d0
                  end if
               else
                  derivative_of_X(n) = (X(i,j,n) - X(i,j-1,n)) * one_over_dy
               end if
            end do

            if (debug2 .and. i .eq. idx .and. j .eq. jdx) then
               print *, "derivative_of_T", derivative_of_T
               print *, "derivative_of_X", derivative_of_X
            end if

c           evaluate diffusion coefficients and other things on the edge
            call CKCPMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), cpms_on_edge)
            call EGSPAR 
     &           (T_on_edge, X_on_edge, Y_on_edge, 
     &           cpms_on_edge, RWRK(egbr), IWRK(egbi))
            call CKMMWY 
     &           (Y_on_edge, IWRK(ckbi), RWRK(ckbr), mmw_on_edge)
c           SUBROUTINE EGSLTDR5 (T, Y, WW, WEG, IWEG, PTC, THETA, D)
            call EGSLTDR5 
     &           (T_on_edge, Y_on_edge, mmw_on_edge, RWRK(egbr), IWRK(egbi),
     &           lambda_prime_on_edge_cgs, rho_theta_on_edge_cgs, rho_Y_D_on_edge_cgs)
            call CKHMS 
     &           (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms_on_edge_cgs)
c           there does not seem to be a chemkin call to get this ratio from T and Y
            pressure_over_rho_on_edge_cgs 
     &           = universal_gas_constant_cgs * T_on_edge / mmw_on_edge

            if (debug2 .and. i .eq. idx .and. j .eq. jdx) then
               print *, "cpms_on_edge", cpms_on_edge
               print *, "mmw_on_edge", mmw_on_edge
               print *, "hms_on_edge_cgs", hms_on_edge_cgs
               print *, "lambda_prime_on_edge_cgs", lambda_prime_on_edge_cgs
               print *, "rho_theta_on_edge_cgs", rho_theta_on_edge_cgs
               print *, "rho_Y_D_on_edge_cgs", rho_Y_D_on_edge_cgs
            end if

c           convert the library values from cgs units to si units

c           per page 98 of the chemkin manual, the enthalpies for the species,
c           hms_on_edge_cgs, have units
c
c                            erg / g
c                 =          (J / 10**7) / (Kg / 10**3)
c                 = 1.0d-4 * J / Kg
c
c           the si unit of thermal conductivity is J / Kg, so the EGLIB values must
c           be multiplied by 1.0d-4 to convert them to si units.
            do n=1,nspecies
               hms_on_edge_si(n) = 1.0d-4 * hms_on_edge_cgs(n)
            end do

c           per page 23 of the EGLIB manual, the partial thermal conductivity, 
c           lambda_prime_on_edge_cgs, has units
c
c                            erg / (cm s K)
c                 =          (J / 10**7) / ((m / 10**2) s K)
c                 = 1.0d-5 * ((W s) / (m s K))
c                 = 1.0d-5 * W / (m K).
c
c           the si unit of thermal conductivity is W / (m K), so the EGLIB value must
c           be multiplied by 1.0d-5 to convert it to si units.
            lambda_prime_on_edge_si = 1.0d-5 * lambda_prime_on_edge_cgs

c           per page 23 in the EGLIB manual, the rescaled thermal diffusion vector,
c           rho_theta_on_edge_cgs, has units
c
c                            g / (cm s)
c                 =          (Kg / 10**3) / ((m / 10**2) s)
c                 = 1.0d-1 * Kg / (m s)
c
c           so the EGLIB values must be multiplied by 1.0d-1 to convert them to si units.
            do n=1,nspecies
               rho_theta_on_edge_si(n) = 1.0d-1 * rho_theta_on_edge_cgs(n)
            end do

c           per page 23 in the EGLIB manual, the rescaled flux diffusion matrix,
c           rho_Y_D_on_edge_cgs, has units
c
c                            g / (cm s)
c                 =          (Kg / 10**3) / ((m / 10**2) s)
c                 = 1.0d-1 * Kg / (m s)
c
c           so the EGLIB values must be multiplied by 1.0d-1 to convert them to si units.
            do ni=1,nspecies
               do nj=1,nspecies
                  rho_Y_D_on_edge_si(ni,nj) = 1.0d-1 * rho_Y_D_on_edge_cgs(ni,nj)
               end do
            end do

c           the cgs unit for pressure is the barye which equals 0.1 Pascal.
c           the si unit for pressure is the Pascal which equals 1 Kg m s**2.

c           the pressure and density ratio, pressure_over_rho_on_edge_cgs, has units
c
c                            barye / (g / cm**3)
c                 =          (1.0d-1 Pa) / ((1.0d-3 Kg) / (1.0d-2 m)**3)
c                 = 1.0d-4 * Pa / (Kg / m**2)
c
c           so the cgs value must be multiplied by 1.0d-4 to convert it to si units.
            pressure_over_rho_on_edge_si = 1.0d-4 * pressure_over_rho_on_edge_cgs 

            if (debug2 .and. i .eq. idx .and. j .eq. jdx) then
               print *, "hms_on_edge_si", hms_on_edge_si
               print *, "lambda_prime_on_edge_si", lambda_prime_on_edge_si
               print *, "rho_theta_on_edge_si", rho_theta_on_edge_si
               print *, "rho_Y_D_on_edge_si", rho_Y_D_on_edge_si
            end if

c           assemble the fluxes
            sum_h_yflux_for_Y = 0
            sum_rho_theta_derivative_of_X = 0
            do n = 1,nspecies
               sum_rho_Y_D_derivative_of_X = 0
               do ell = 1,nspecies
                  sum_rho_Y_D_derivative_of_X 
     &                 = sum_rho_Y_D_derivative_of_X 
     &                 + rho_Y_D_on_edge_si(n,ell) 
     &                 * derivative_of_X(ell)
               end do
               yflux_for_Y(i,j,n) 
     &              = - area_for_y
     &              * (sum_rho_Y_D_derivative_of_X
     &              + Y_on_edge(n) * rho_theta_on_edge_si(n) * derivative_of_T / T_on_edge)
               sum_h_yflux_for_Y 
     &              = sum_h_yflux_for_Y 
     &              + hms_on_edge_si(n) 
     &              * yflux_for_Y(i,j,n)
               sum_rho_theta_derivative_of_X 
     &              = sum_rho_theta_derivative_of_X 
     &              + rho_theta_on_edge_si(n) 
     &              * derivative_of_X(n)
            end do
            yflux_for_H(i,j)
     &           = sum_h_yflux_for_Y
     &           - area_for_y
     &           * (lambda_prime_on_edge_si * derivative_of_T
     &           + pressure_over_rho_on_edge_si * sum_rho_theta_derivative_of_X)
            if (debug2 .and. i .eq. idx .and. j .eq. jdx) then
               print *, "sum_h_yflux_for_Y", sum_h_yflux_for_Y
               print *, "- area_for_y", - area_for_y
               print *, "lambda_prime_on_edge_si", lambda_prime_on_edge_si
               print *, "derivative_of_T", derivative_of_T
               print *, "pressure_over_rho_on_edge_si", pressure_over_rho_on_edge_si
               print *, "sum_rho_theta_derivative_of_X", sum_rho_theta_derivative_of_X
               print *, "term 2", - area_for_y * lambda_prime_on_edge_si * derivative_of_T
               print *, "term 3", - area_for_y * pressure_over_rho_on_edge_si * sum_rho_theta_derivative_of_X
               print *, "yflux_for_H(i,j)", yflux_for_H(i,j)
            end if

            end do
         end do

      if (debug) call print_patch_n ("yflux_for_H", yflux_for_H, DIMS(yflux_for_H), 1, idx, jdx)
      if (debug) call print_patch_n ("yflux_for_Y", yflux_for_Y, DIMS(yflux_for_Y), nspecies, idx, jdx)

c///
c///  form the updates weighted by scale_factor / vol
c///

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            div_of_flux_for_H(i, j) 
     &           = scale_factor_over_volume
     &           * ((xflux_for_H(i+1,j) - xflux_for_H(i,j)) 
     &           + (yflux_for_H(i,j+1) - yflux_for_H(i,j)))
            do n=1,nspecies
               div_of_flux_for_Y(i, j, n) 
     &              = scale_factor_over_volume
     &              * ((xflux_for_Y(i+1,j,n) - xflux_for_Y(i,j,n))
     &              + (yflux_for_Y(i,j+1,n) - yflux_for_Y(i,j,n)))
            end do
         end do
      end do

      if (debug) call print_patch_n ("div_of_flux_for_H", div_of_flux_for_H, DIMS(div_of_flux_for_H), 1, idx, jdx)
      if (debug) call print_patch_n ("div_of_flux_for_Y", div_of_flux_for_Y, DIMS(div_of_flux_for_Y), nspecies, idx, jdx)

c///
c///  release temporary arrays (not strictly necessary)
c///

      deallocate (contiguous_X)
      deallocate (contiguous_Y)
      deallocate (cpms_on_edge)
      deallocate (derivative_of_X)
      deallocate (hms_on_edge_cgs)
      deallocate (hms_on_edge_si)
      deallocate (residual)
      deallocate (rho_theta_on_edge_cgs)
      deallocate (rho_theta_on_edge_si)
      deallocate (rho_Y_D_on_edge_cgs)
      deallocate (rho_Y_D_on_edge_si)
      deallocate (T)
      deallocate (X)
      deallocate (X_on_edge)
      deallocate (Y)
      deallocate (Y_on_edge)

c///
c///  end of FORT_RK_MULTICOMPONENT
c///

      if (debug) print *, "JFG: at bottom of FORT_RK_MULTICOMPONENT"
c      if (debug) stop

      return
      end

      subroutine print_patch_n (
     &     title,
     &     array, lo1, lo2, hi1, hi2, ncomps,
     &     i_loc, j_loc)

      character title*(*)
      integer i
      integer i_loc
      integer j
      integer j_loc
      integer lo1, lo2, hi1, hi2
      integer n
      integer ncomps
      REAL_T array(lo1 : hi1, lo2 : hi2, ncomps)

      write (6, 60000) title
      if (lo1 .le. i_loc-1 .and. i_loc+1 .le. hi1 .and.
     &    lo2 .le. j_loc-1 .and. j_loc+1 .le. hi2) then
         do n=1,ncomps
            write (6, 60001) (i, i=i_loc-1,i_loc+1)
            write (6, 60002) (j, (array(i, j, n), i=i_loc-1,i_loc+1), j=j_loc+1,j_loc-1,-1)
         end do
         write (6, "()")
      else
         write (6, "()")
         write (6, 60000) 
     &        "the print region is not in the range of this box"
         write (6, "()")
      end if

60000 format (4x, a)
60001 format (/6x, 3i15)
60002 format (1p, 1x, i5, 3d15.6 /1x, i5, 3d15.6 /1x, i5, 3d15.6)

      return
      end
