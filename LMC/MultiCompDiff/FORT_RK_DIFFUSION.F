c#undef BL_LANG_CC
c#ifndef BL_LANG_FORT
c#define BL_LANG_FORT
c#endif

c#define BL_SPACEDIM 2

c#include "ArrayLim.H"

c     to do:
c     must adjust diffusion coefficients by thickness factors
c     unsure of LMC support for r-z coordinates

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

      subroutine FORT_RK_DIFFUSION
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
      dimension X (:)
      dimension Y (:,:,:)

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
      allocate (X (nspecies))
      allocate (Y (DIMV(state), nspecies))

c     offset of species in the state
      offset = index_of_firstY - 1

c     copy temperature from the state
      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            temperature(i, j) = state(i, j, index_of_T)
         end do
      end do

c     calculate H and Y from rho, rho H and rho Y in the state
      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            one_over_rho = 1.0d0 / state(i, j, index_of_rho)
            H(i, j) = one_over_rho * state(i, j, index_of_rhoH)
            do n=1,nspecies
               rho_Y(i, j, n) = state(i, j, offset + n)
               Y(i, j, n) = one_over_rho * state(i, j, offset + n)
            end do
         end do
      end do

c     get ambient pressure in cgs units
c     one_atmosphere_cgs = 1013250.0 (dyne / cm**2)
c     one_atmosphere_si = 101325.0 (pascal)
      call CKRP (IWRK(ckbi), RWRK(ckbr), RU, RUC, one_atmosphere_cgs)
      call FORT_GETPAMB (ambient_pressure_si, dpdt_factor)
      one_atmosphere_si = FORT_P1ATMMKS ()
      ambient_pressure_cgs = one_atmosphere_cgs * ambient_pressure_si / one_atmosphere_si

c     calculate temperature, rho_D, and lambda

      call CKWT (IWRK(ckbi),RWRK(ckbr),one_over_wt)
      do n=1,nspecies
         one_over_wt(n) = 1.0d0 / one_over_wt(n)
      end do

      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            do n=1,nspecies
               contiguous_Y(n) = Y(i,j,n)
            end do

c           calculate temperature consistent with H and Y
            iterations = FORT_TfromHYpt 
     &           (temperature(i, j), 
     &           H(i, j), 
     &           contiguous_Y,
     &           maximum_error,
     &           maximum_iterations,
     &           residual)

            floored_T = MAX (temperature(i,j), TMIN_TRANS) 

            CALL CKCPMS (floored_T, IWRK(ckbi), RWRK(ckbr), cpms)
            CALL CKYTX (contiguous_Y, IWRK(ckbi), RWRK(ckbr), X)
            CALL EGSPAR 
     &           (floored_T, X, contiguous_Y, cpms, RWRK(egbr), 
     &           IWRK(egbi))

            CALL CKMMWY (contiguous_Y, IWRK(ckbi), RWRK(ckbr), mmw)
            CALL EGSV1 
     &           (ambient_pressure_cgs, floored_T, contiguous_Y, mmw, RWRK(egbr), 
     &           Upsilon)

c           note, this is the practice in the FORT_SPECTEMPVISC routine that underlies 
c           calcDiffusivity, but since rho is in the state, why recompute it? 
            CALL CKRHOY (ambient_pressure_cgs, floored_T, contiguous_Y, IWRK(ckbi),
     &           RWRK(ckbr), rho_cgs)
            do n=1,nspecies
c              per page 17 in the EGLIB manual, Upsilon has units (cm**2 / s), so
c              rho Upsilon has units (g / cm s).  per equation 12 in the EGLIB manual, 
c              the mixture averaged diffusion coefficient is rho Upsilon times the 
c              dimensionless mmw / wt(n).  thus rho times the diffusion coefficient 
c              with units (g / (cm s)) = 10 (Kg / (m s)) must be divided by 10 to obtain
c              si units.
               rho_D(i,j,n) = rho_cgs * mmw * one_over_wt(n) * Upsilon(n) * 0.1d0
            end do

c           the following averaging is equation 50 in the chemkin transport package 
c           manual.  the EGLIB manual does not mention the averaging, but it does have
c           provisions to compute it through the alpha parameter in its equation 17.
            alpha = 1.0d0
            CALL EGSL1 (alpha, floored_T, X, RWRK(egbr), plus_value)
            alpha = - 1.0d0
            CALL EGSL1 (alpha, floored_T, X, RWRK(egbr), minus_value)
c           per page 15 of the EGLIB manual, these thermal conductivities are returned
c           with units (erg / (cm s K)) = ((J / 10**7) / ((m / 10**2) s K)) =
c           = ((W s) / (10**5 m s K)) = (W / (m K 10**5)).  the si unit of thermal
c           conductivity is (W / (M K)), so these values must be divided by 10**5
c           to convert them to si units.   multiplication by 0.5 averages the values.
            lambda(i,j) = 0.5d-5 * (plus_value + minus_value)
         end do
      end do

c     fluxes in x direction are from lo(1) to hi(1)+1
      half_over_dx = 0.5d0 / dx(1)
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
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
     &           * (temperature(i-1,j) 
     &            + temperature(i,j))
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
     &           * (temperature(i,j) - temperature(i-1,j))
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
         one_over_dx = 1.0d0 / dx(1)
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
               temperature(i-1,j) = state(i-1,j,index_of_T)
               T_on_edge = temperature(i-1,j)
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
     &                 * (- (8.0d0 / 3.0d0) * temperature(i-1,j)
     &                 + 3 * temperature(i,j)
     &                 - (1.0d0 / 3.0d0) * temperature(i+1,j))
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
               temperature(i+1,j) = state(i+1,j,index_of_T)
               T_on_edge = temperature(i+1,j)
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
     &              * ((1.0d0 / 3.0d0) * temperature(i-1,j)
     &              - 3 * temperature(i,j)
     &              + (8.0d0 / 3.0d0) * temperature(i+1,j))
            else if (.not. internal) then
               xflux_for_H(i+1,j) = 0
            end if
         end do
      end if

c     fluxes in y direction are from lo(2) to hi(2)+1
      half_over_dy = 0.5d0 / dx(2)
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
     &           * (temperature(i,j-1) 
     &            + temperature(i,j))
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
     &           * (temperature(i,j) - temperature(i,j-1))
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
         one_over_dy = 1.0d0 / dx(2)
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
               temperature(i,j-1) = state(i,j-1,index_of_T)
               T_on_edge = temperature(i,j-1)
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
     &                 * (- (8.0d0 / 3.0d0) * temperature(i,j-1)
     &                 + 3 * temperature(i,j)
     &                 - (1.0d0 / 3.0d0) * temperature(i,j+1))
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
               temperature(i,j+1) = state(i,j+1,index_of_T)
               T_on_edge = temperature(i,j+1)
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
     &              * ((1.0d0 / 3.0d0) * temperature(i,j-1)
     &              - 3 * temperature(i,j)
     &              + (8.0d0 / 3.0d0) * temperature(i,j+1))
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
