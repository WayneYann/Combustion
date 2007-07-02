c#undef BL_LANG_CC
c#ifndef BL_LANG_FORT
c#define BL_LANG_FORT
c#endif

c#define BL_SPACEDIM 2

c#include "ArrayLim.H"

c     to do:
c     must adjust diffusion coefficients by thickness factors
c     unsure of the signs and scaling wanted for fluxes and updates:
c        at present the updates (but not the fluxes) are weighted by 1/vol
c        so they give the change to the pointwise value
c     unsure of LMC support for r-z coordinates

c     arguments are alphabetical, mostly:
c
c     domain_lo, domain_hi,             ! INPUT limits of valid region of the domain
c     lo, hi,                           ! INPUT limits of valid region of the box
c     areax, DIMS(areax),               ! INPUT areas of the faces perpendicular to x axis
c     areay, DIMS(areay),               ! INPUT areas of the faces perpendicular to y axis
c     bc,                               ! INPUT boundary condition array for all comps
c     dx,                               ! INPUT physical dimensions of grid cells
c     index_of_firstY,                  ! INPUT index of rho Y for the first species in the state
c     index_of_lastY,                   ! INPUT index of rho Y for the last species in the state
c     index_of_rho,                     ! INPUT index of rho in the state
c     index_of_rhoH,                    ! INPUT index of rho H in the state
c     index_of_T                        ! INPUT index of T in the state
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
     &     dx,
     &     index_of_firstY,
     &     index_of_lastY,
     &     index_of_rho,
     &     index_of_rhoH,
     &     index_of_T,
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
c#include "CONSTANTS.H"
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
      REAL_T areax(DIMV(areax))
      REAL_T areay(DIMV(areay))
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
      integer i
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
      REAL_T contiguous_Y
      REAL_T correction
      REAL_T cpms
      REAL_T floored_T
      REAL_T half_over_dx
      REAL_T half_over_dy
      REAL_T hms 
      REAL_T hms_scale_factor
      REAL_T lambda
      REAL_T minus_value
      REAL_T mmw
      REAL_T neg_one_over_volume
      REAL_T one_over_dx
      REAL_T one_over_rho
      REAL_T one_over_wt
      REAL_T Patm
      REAL_T P1ATM
      REAL_T plus_value
      REAL_T pressure
      REAL_T rho
      REAL_T rho_D
      REAL_T RU
      REAL_T RUC
      REAL_T sum
      REAL_T T_on_edge
      REAL_T Upsilon
      REAL_T X 
      REAL_T Y 
      allocatable contiguous_Y
      allocatable cpms
      allocatable hms
      allocatable lambda
      allocatable one_over_wt
      allocatable rho_D
      allocatable Upsilon
      allocatable X
      allocatable Y
      dimension contiguous_Y (:)
      dimension cpms (:)
      dimension hms (:)
      dimension lambda (:,:)
      dimension one_over_wt (:)
      dimension rho_D (:,:,:)
      dimension Upsilon (:)
      dimension X (:)
      dimension Y (:,:,:)

c     allocate all variably dimensioned local arrays
      allocate (contiguous_Y (nspecies))
      allocate (cpms (nspecies))
      allocate (hms (nspecies))
      allocate (lambda (DIMV(state)))
      allocate (one_over_wt (nspecies))
      allocate (rho_D (DIMV(state), nspecies))
      allocate (Upsilon (nspecies))
      allocate (X (nspecies))
      allocate (Y (DIMV(state), nspecies))

c     offset of species in the state
      offset = index_of_firstY - 1

c     calculate Y from rho Y in the state
      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            one_over_rho = 1.0d0 / state(i, j, index_of_rho)
            do n=1,nspecies
               Y(i, j, n) = one_over_rho * state(i, j, offset + n)
            end do
         end do
      end do

c     calculate rho_D and lambda

      CALL CKRP (IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      pressure = Patm * P1ATM

      call CKWT (IWRK(ckbi),RWRK(ckbr),one_over_wt)
      do n=1,nspecies
         one_over_wt(n) = 1.0d0 / one_over_wt(n)
      end do

      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            do n=1,nspecies
               contiguous_Y(n) = Y(i,j,n)
            end do
            floored_T = MAX (state(i,j,index_of_T), TMIN_TRANS) 

            CALL CKCPMS (floored_T, IWRK(ckbi), RWRK(ckbr), cpms)
            CALL CKYTX (contiguous_Y, IWRK(ckbi), RWRK(ckbr), X)
            CALL EGSPAR 
     &           (floored_T, X, contiguous_Y, cpms, RWRK(egbr), 
     &           IWRK(egbi))

            CALL CKMMWY (contiguous_Y, IWRK(ckbi), RWRK(ckbr), mmw)
            CALL EGSV1 
     &           (pressure, floored_T, contiguous_Y, mmw, RWRK(egbr), 
     &           Upsilon)

c           since rho is in the state, why recompute it? 
c           note, this is the practice in the FORT_SPECTEMPVISC routine
c           that underlies calcDiffusivity.
            CALL CKRHOY (pressure, floored_T, contiguous_Y, IWRK(ckbi),
     &           RWRK(ckbr), rho)
            do n=1,nspecies
c              per equation 12 in the EGLIB manual.
c              multiplication by 0.1d0 converts from cgs to mks.
               rho_D(i,j,n) 
     &              = rho * mmw * one_over_wt(n) * Upsilon(n) * 0.1d0
            end do

c           this averaging is equation 50 in the chemkin transport package manual.
c           note, the EGLIB manual does not mention the averaging, but it does have
c           provisions to compute it through the alpha parameter in equation 17.
            alpha = 1.0d0
            CALL EGSL1 (alpha, floored_T, X, RWRK(egbr), plus_value)
            alpha = - 1.0d0
            CALL EGSL1 (alpha, floored_T, X, RWRK(egbr), minus_value)
c           multiplication by 1.0e-5 converts from cgs to mks
c           multiplication by 0.5 averages the values
            lambda(i,j) = 0.5d-5 * (plus_value + minus_value)
         end do
      end do

c     multiplication by hms_scale_factor converts from cgs to MKS (1 erg/g = 1.e-4 J/kg)
      hms_scale_factor = 1.0d-4

c     fluxes in x direction
      half_over_dx = 0.5d0 / dx(1)
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
c           species fluxes
            sum = 0
            do n=1,nspecies
               xflux_for_Y(i,j,n) 
     &              = half_over_dx 
     &              * areax(i,j) 
     &              * (rho_D(i-1,j,n) + rho_D(i,j,n)) 
     &              * (Y(i,j,n) - Y(i-1,j,n))
               sum = sum + xflux_for_Y(i,j,n)
            end do
c           flux correction
            correction = - sum / nspecies
            do n=1,nspecies
               xflux_for_Y(i,j,n) = xflux_for_Y(i,j,n) + correction
            end do
c           enthalpy flux
            T_on_edge = 0.5d0
     &           * (state(i-1,j,index_of_T) 
     &            + state(i,j,index_of_T))
            call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
            sum = 0
            do n=1,nspecies
               sum = sum + xflux_for_Y(i,j,n) * hms(n) * hms_scale_factor
            end do
            xflux_for_H(i,j) 
     &           = sum
     &           + half_over_dx
     &           * areax(i,j) 
     &           * (lambda(i-1,j) + lambda(i,j)) 
     &           * (state(i,j,index_of_T) - state(i-1,j,index_of_T))
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
               sum = 0
               do n=1,nspecies
                  xflux_for_Y(i,j,n) 
     &                 = one_over_dx 
     &                 * areax(i,j) 
     &                 * rho_D(i-1,j,n)
     &                 * (- (8.0d0 / 3.0d0) * Y(i-1,j,n)
     &                 + 3 * Y(i,j,n)
     &                 - (1.0d0 / 3.0d0) * Y(i+1,j,n))
                  sum = sum + xflux_for_Y(i,j,n) 
               end do
c              flux correction
               correction = - sum / nspecies
               do n=1,nspecies
                  xflux_for_Y(i,j,n) = xflux_for_Y(i,j,n) + correction
               end do
            else if (.not. internal) then
               do n=1,nspecies
                  xflux_for_Y(i,j,n) = 0
               end do
            end if
            if (external_for_T) then          
               T_on_edge = state(i-1,j,index_of_T)
               call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
               sum = 0
               do n=1,nspecies
                  sum = sum + xflux_for_Y(i,j,n) * hms(n) * hms_scale_factor
               end do
               xflux_for_H(i,j) 
     &              = sum
     &              + one_over_dx
     &              * areax(i,j) 
     &              * lambda(i-1,j)
     &                 * (- (8.0d0 / 3.0d0) * state(i-1,j,index_of_T)
     &                 + 3 * state(i,j,index_of_T)
     &                 - (1.0d0 / 3.0d0) * state(i+1,j,index_of_T))
            else if (.not. internal) then
               xflux_for_H(i,j) = 0
            end if
         end do

c     redo the fluxes at high x domain boundary
      else if (domain_hi(1) .eq. hi(1)) then
c        check the bc flags
         external_for_T = bc(1, 2, index_of_T) .eq. EXT_DIR
         internal = bc(1, 2, index_of_T) .eq. INT_DIR
         external_for_Y = .true.
         external_check = .true.
         internal_check = .true.
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
               sum = 0
               do n=1,nspecies
                  xflux_for_Y(i+1,j,n) 
     &                 = one_over_dx 
     &                 * areax(i+1,j) 
     &                 * rho_D(i+1,j,n)
     &                 * ((1.0d0 / 3.0d0) * Y(i-1,j,n)
     &                 - 3 * Y(i,j,n)
     &                 + (8.0d0 / 3.0d0) * Y(i+1,j,n))
                  sum = sum + xflux_for_Y(i+1,j,n) 
               end do
c              flux correction
               correction = - sum / nspecies
               do n=1,nspecies
                  xflux_for_Y(i+1,j,n) = xflux_for_Y(i+1,j,n) + correction
               end do
            else if (.not. internal) then
               do n=1,nspecies
                  xflux_for_Y(i+1,j,n) = 0
               end do
            end if
            if (external_for_T) then          
               T_on_edge = state(i+1,j,index_of_T)
               call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
               sum = 0
               do n=1,nspecies
                  sum = sum + xflux_for_Y(i+1,j,n) * hms(n) * hms_scale_factor
               end do
               xflux_for_H(i+1,j) 
     &              = sum
     &              + one_over_dx
     &              * areax(i+1,j) 
     &              * lambda(i+1,j)
     &              * ((1.0d0 / 3.0d0) * state(i-1,j,index_of_T)
     &              - 3 * state(i,j,index_of_T)
     &              + (8.0d0 / 3.0d0) * state(i+1,j,index_of_T))
            else if (.not. internal) then
               xflux_for_H(i+1,j) = 0
            end if
         end do
      end if

c     fluxes in y direction
      half_over_dx = 0.5d0 / dx(2)
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
c           species fluxes
            sum = 0
            do n=1,nspecies
               yflux_for_Y(i,j,n) 
     &              = half_over_dy
     &              * areay(i,j) 
     &              * (rho_D(i,j-1,n) + rho_D(i,j,n)) 
     &              * (Y(i,j,n) - Y(i,j-1,n))
               sum = sum + yflux_for_Y(i,j,n)
            end do
c           flux correction
            correction = - sum / nspecies
            do n=1,nspecies
               yflux_for_Y(i,j,n) = yflux_for_Y(i,j,n) + correction
            end do
c           enthalpy flux
            T_on_edge = 0.5d0
     &           * (state(i,j-1,index_of_T) 
     &            + state(i,j,index_of_T))
            call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
            sum = 0
            do n=1,nspecies
               sum = sum + xflux_for_Y(i,j,n) * hms(n) * hms_scale_factor
            end do
            yflux_for_H(i,j) 
     &           = sum
     &           + half_over_dx
     &           * areax(i,j) 
     &           * (lambda(i,j-1) + lambda(i,j)) 
     &           * (state(i,j,index_of_T) - state(i,j-1,index_of_T))
         end do
      end do

c     redo the fluxes at low y domain boundary
c     redo the fluxes in y direction at domain boundaries
      if (domain_lo(2) .eq. lo(2)) then
         j = lo(2)
         do i=lo(1),hi(1)
c           do species first because their fluxes appear in the enthalpy case
            do n=1,nspecies
               if (bc(2, 1, offset + n) .eq. EXT_DIR) then
c                 unsymmetric stencil for dirichlet conditions
                  yflux_for_Y(i,j,n) 
     &                 = half_over_dx 
     &                 * areax(i,j) 
     &                 * (rho_D(i,j-1,n) + rho_D(i,j,n)) 
     &                 * (- (8.0d0 / 3.0d0) * Y(i,j-1,n)
     &                 + 3 * Y(i,j,n)
     &                 - (1.0d0 / 3.0d0) * Y(i,j+1,n))
               else if (bc(2, 1, offset + n) .ne. INT_DIR) then
c                 no flux in all other cases (INT_DIR should never occur 
c                 at a domain boundary but it does not hurt to check)
                  yflux_for_Y(i,j,n) = 0
               end if
            end do
c           now the enthalpy case
            if (bc(2, 1, index_of_T) .eq. EXT_DIR) then
c              unsymmetric stencil for dirichlet conditions
               sum = half_over_dx
     &              * areay(i,j) 
     &              * (lambda(i,j-1) + lambda(i,j)) 
     &              * (- (8.0d0 / 3.0d0) * state(i,j-1,index_of_T)
     &              + 3 * state(i,j,index_of_T)
     &              - (1.0d0 / 3.0d0) * state(i,j+1,index_of_T))
               T_on_edge = 0.5d0
     &              * (state(i,j-1,index_of_T) 
     &              + state(i,j,index_of_T))
               call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
               do n=1,nspecies
                  sum = sum
     &                 + yflux_for_Y(i,j,n)
     &                 * hms(n)
     &                 * hms_scale_factor
               end do
               yflux_for_H(i,j) = sum
            else if (bc(2, 1, index_of_T) .ne. INT_DIR) then
c              no flux in all other cases (INT_DIR should never occur 
c              at a domain boundary but it does not hurt to check).
c              note the code logic assumes if neither EXT_DIR nor INT_DIR 
c              for T, then the same is true for all the species.
               yflux_for_H(i,j) = 0
            end if
         end do
      else if (domain_hi(2) .eq. hi(2)) then
         j = hi(2)
         do i=lo(1),hi(1)
c           do species first because their fluxes appear in the enthalpy case
            do n=1,nspecies
               if (bc(2, 2, offset + n) .eq. EXT_DIR) then
c                 unsymmetric stencil for dirichlet conditions
                  yflux_for_Y(i,j+1,n) 
     &                 = half_over_dx 
     &                 * areax(i,j) 
     &                 * (rho_D(i,j-1,n) + rho_D(i,j,n)) 
     &                 * ((1.0d0 / 3.0d0) * Y(i,j-1,n)
     &                 - 3 * Y(i,j,n)
     &                 + (8.0d0 / 3.0d0) * Y(i,j+1,n))
               else if (bc(2, 2, offset + n) .ne. INT_DIR) then
c                 no flux in all other cases (INT_DIR should never occur 
c                 at a domain boundary but it does not hurt to check)
                  yflux_for_Y(i,j+1,n) = 0
               end if
            end do
c           now the enthalpy case
            if (bc(2, 2, index_of_T) .eq. EXT_DIR) then
c              unsymmetric stencil for dirichlet conditions
               sum = half_over_dx
     &              * areay(i,j) 
     &              * (lambda(i,j-1) + lambda(i,j)) 
     &              * ((1.0d0 / 3.0d0) * state(i,j-1,index_of_T)
     &              - 3 * state(i,j,index_of_T)
     &              + (8.0d0 / 3.0d0) * state(i,j+1,index_of_T))
               T_on_edge = half 
     &              * (state(i,j-1,index_of_T) 
     &              + state(i,j,index_of_T))
               call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
               do n=1,nspecies
                  sum = sum
     &                 + yflux_for_Y(i,j+1,n)
     &                 * hms(n)
     &                 * hms_scale_factor
               end do
               yflux_for_H(i,j+1) = sum
            else if (bc(2, 2, index_of_T) .ne. INT_DIR) then
c              no flux in all other cases (INT_DIR should never occur 
c              at a domain boundary but it does not hurt to check).
c              note the code logic assumes if neither EXT_DIR nor INT_DIR 
c              for T, then the same is true for all the species.
               yflux_for_H(i,j+1) = 0
            end if
         end do
      end if

c     volume weighted negative sum of the fluxes
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            neg_one_over_volume = - 1.0d0 / volume(i,j)
            update_for_H(i, j) 
     &           = neg_one_over_volume
     &           * (xflux_for_H(i+1,j) - xflux_for_H(i,j)
     &           + yflux_for_H(i,j+1) - yflux_for_H(i,j))
            do n=1,nspecies
               update_for_Y(i, j, n) 
     &              = neg_one_over_volume
     &              * (xflux_for_Y(i+1,j,n) - xflux_for_Y(i,j,n)
     &              + yflux_for_Y(i,j+1,n) - yflux_for_Y(i,j,n))
            end do
         end do
      end do

      deallocate (contiguous_Y)
      deallocate (cpms)
      deallocate (hms)
      deallocate (lambda)
      deallocate (one_over_wt)
      deallocate (rho_D)
      deallocate (Upsilon)
      deallocate (X)
      deallocate (Y)

      return
      end
