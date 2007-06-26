#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define BL_SPACEDIM 2

#include "ArrayLim.H"

      subroutine FORT_RK_DIFFUSION
     &     (
     &     lo, hi,               ! INPUT limits of valid region
     &     areax, DIMS(areax),   ! INPUT areas of the faces perpendicular to x axis
     &     areay, DIMS(areay),   ! INPUT areas of the faces perpendicular to y axis
     &     dx,                   ! INPUT physical dimensions of grid cells
     &     state, DIMS(state),   ! INPUT all variables in the state
     &     update_for_H, DIMS(update_for_H), ! OUTPUT divergences of the fluxes
     &     update_for_Y, DIMS(update_for_Y), ! OUTPUT divergences of the fluxes
     &     xflux_for_H, DIMS(xflux_for_H),   ! OUTPUT x fluxes for enthalpy
     &     yflux_for_H, DIMS(yflux_for_H),   ! OUTPUT y fluxes for enthalpy
     &     xflux_for_Y, DIMS(xflux_for_Y),   ! OUTPUT x fluxes for species
     &     yflux_for_Y, DIMS(yflux_for_Y),   ! OUTPUT y fluxes for species
     &     vol, DIMS(vol),       ! INPUT volumes of the cells
     &     ncomps,               ! INPUT total number of components in the state
     &     nspecies,             ! INPUT total number of species in the state
     &     index_of_firstY,      ! INPUT index of rho Y for the first species in the state
     &     index_of_lastY,       ! INPUT index of rho Y for the last species in the state
     &     index_of_rho,         ! INPUT index of rho in the state
     &     index_of_rhoH,        ! INPUT index of rho H in the state
     &     index_of_T            ! INPUT index of T in the state
     &     )

      implicit none

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "HEATTRANSFER_F.H"
#include "ChemDriver_F.H"
#include "cdwrk.H"

c     to do
c     must adjust diffusion coefficients by thickness factors
c     must include fixes for boundary conditions

c     arguments
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer DIMDEC(areax)
      integer DIMDEC(areay)
      integer DIMDEC(state)
      integer DIMDEC(update_for_H)
      integer DIMDEC(update_for_Y)
      integer DIMDEC(xflux_for_H)
      integer DIMDEC(yflux_for_H)
      integer DIMDEC(xflux_for_Y)
      integer DIMDEC(yflux_for_Y)
      integer DIMDEC(vol)
      integer ncomps
      integer nspecies
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
      REAL_T vol(DIMV(vol))

c     local variables
      integer i, j, n, offset
      REAL_T alpha
      REAL_T contiguous_Y
      REAL_T cpms
      REAL_T floored_T
      REAL_T half_over_dx
      REAL_T half_over_dy
      REAL_T half_over_hundred_thousand
      REAL_T hms ! species enthalpies at an edge
      REAL_T lambda
      REAL_T minus_value
      REAL_T mmw
      REAL_T neg_one_over_vol
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
      REAL_T scale_factor
      REAL_T sum
      REAL_T T_on_edge
      REAL_T Upsilon
      REAL_T X ! mole fractions at cell centers
      REAL_T Y ! mass fractions at cell centers
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

c     calculate Y from rho Y in the state
      allocate (Y (DIMV(state), nspecies))
      offset = index_of_firstY - 1
      do j=ARG_L2(state),ARG_H2(state)
         do i=ARG_L1(state),ARG_H1(state)
            one_over_rho = one / state(j, i, index_of_rho)
            do n=1,nspecies
               Y(j, i, n) = one_over_rho * state(j, i, offset + n)
            end do
         end do
      end do

c     calculate rho_D and lambda
      allocate (contiguous_Y (nspecies))
      allocate (lambda (DIMV(state)))
      allocate (one_over_wt (nspecies))
      allocate (rho_D (DIMV(state), nspecies))
      allocate (Upsilon (nspecies))
      allocate (X (nspecies))

      CALL CKRP (IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      pressure = Patm * P1ATM

      call CKWT (IWRK(ckbi),RWRK(ckbr),one_over_wt)
      do n=1,nspecies
         one_over_wt(n) = one / one_over_wt(n)
      end do

      half_over_hundred_thousand = half / 100000
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
c           that does much of the work for calcDiffusivity.
            CALL CKRHOY (pressure, floored_T, contiguous_Y, IWRK(ckbi),
     &           RWRK(ckbr), rho)
            do n=1,nspecies
c              per equation 12 in the EGLIB manual
c              multiplication by tenth converts from cgs to mks
               rho_D(i,j,n) 
     &              = rho * mmw * one_over_wt(n) * Upsilon(n) * tenth
            end do

c           this averaging per equation 50 in the chemkin transport package manual.
c           note, the EGLIB manual does not mention the averaging, but it does have
c           provisions to compute it through the alpha parameter in equation 17
            alpha = one
            CALL EGSL1 (alpha, floored_T, X, RWRK(egbr), plus_value)
            alpha = - one
            CALL EGSL1 (alpha, floored_T, X, RWRK(egbr), minus_value)
c           multiplication by 1.0e-5 converts from cgs to mks
            lambda(i,j) 
     &           = half_over_hundred_thousand
     &           * (plus_value + minus_value)
         end do
      end do

      allocate (hms (nspecies))
c     multiplication by scale_factor converts from cgs to MKS (1 erg/g = 1.e-4 J/kg)
      scale_factor = tenth**4

c     fluxes in x direction
      half_over_dx = half / dx(1)
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            sum = half_over_dy
     &           * areay(i,j) 
     &           * (lambda(i-1,j) + lambda(i,j)) 
     &           * (state(i,j,index_of_T) - state(i-1,j,index_of_T))
            T_on_edge = half 
     &           * (state(i-1,j,index_of_T) 
     &            + state(i,j,index_of_T))
            call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
            do n=1,nspecies
               xflux_for_Y(i,j,n) 
     &              = half_over_dx 
     &              * areax(i,j) 
     &              * (rho_D(i-1,j,n) + rho_D(i,j,n)) 
     &              * (Y(i,j,n) - Y(i-1,j,n))
c     per Marc, in the enthalpy flux, we want to preserve the mass
c     flux as a subexpression, so we simply multiply by an egde h.
               sum = sum
     &              + xflux_for_Y(i,j,n)
     &              * hms(n)
     &              * scale_factor
            end do
            xflux_for_H(i,j) = sum
         end do
      end do

c     fluxes in y direction
      half_over_dy = half / dx(2)
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            sum = half_over_dy
     &           * areay(i,j) 
     &           * (lambda(i,j-1) + lambda(i,j)) 
     &           * (state(i,j,index_of_T) - state(i,j-1,index_of_T))
            T_on_edge = half 
     &           * (state(i,j-1,index_of_T) 
     &            + state(i,j,index_of_T))
            call CKHMS (T_on_edge, IWRK(ckbi), RWRK(ckbr), hms)
            do n=1,nspecies
               yflux_for_Y(i,j,n) 
     &              = half_over_dy
     &              * areay(i,j) 
     &              * (rho_D(i,j-1,n) + rho_D(i,j,n)) 
     &              * (Y(i,j,n) - Y(i,j-1,n))
c     per Marc, in the enthalpy flux, we want to preserve the mass
c     flux as a subexpression, so we simply multiply by an egde h.
               sum = sum
     &              + yflux_for_Y(i,j,n)
     &              * hms(n)
     &              * scale_factor
            end do
            yflux_for_H(i,j) = sum
         end do
      end do

c     volume weighted negative sum of the fluxes
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            neg_one_over_vol = - one / vol(i,j)
            update_for_H(i, j) 
     &           = neg_one_over_vol
     &           * (xflux_for_H(i+1,j) - xflux_for_H(i,j)
     &           + yflux_for_H(i,j+1) - yflux_for_H(i,j))
            do n=1,nspecies
               update_for_Y(i, j, n) 
     &              = neg_one_over_vol
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
