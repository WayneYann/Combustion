c#undef BL_LANG_CC
c#ifndef BL_LANG_FORT
c#define BL_LANG_FORT
c#endif

c#define BL_SPACEDIM 2

c#include "ArrayLim.H"

c     estimate the time step for the Runge-Kutta implementation of the diffusion
c     update.  this estimation is based on the diffusion coefficients of the 
c     species into the mixture, independent of whether or not the "operator" is 
c     based on the mixture averaged or the multicomponent formulation with Soret 
c     and Dufour effects.
c
c                               space_step ** 2
c                 time_step  =  ---------------
c                               D_mix * 2 * DIM
c
c     the minimum of this formula is taken over all cells and species

c     to do:
c     1) must adjust diffusion coefficients by thickness factors
c     2) unsure of LMC support for r-z coordinates
c     3) could include an estimate for thermal diffusion

c     arguments are alphabetical, mostly:
c
c     lo, hi             ! INPUT limits of valid region of the box
c     dx                 ! INPUT physical dimensions of grid cells
c     index_of_firstY    ! INPUT index of rho Y for the first species in the state
c     index_of_lastY     ! INPUT index of rho Y for the last species in the state
c     index_of_rho       ! INPUT index of rho in the state
c     index_of_rhoH      ! INPUT index of rho H in the state
c     index_of_T         ! INPUT index of T in the state
c     maximum_error      ! INPUT maximum error in calculation of T
c     maximum_iterations ! INPUT maximum iterations in calculation of T
c     ncomps             ! INPUT total number of components in the state
c     nspecies           ! INPUT total number of species in the state
c     smallest_time_step ! OUTPUT smallest time step over all cells and species
c     state, DIMS(state) ! INPUT all variables in the state

      subroutine FORT_RK_STEP_SELECTION
     &     (
     &     lo, hi,
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
     &     smallest_time_step,
     &     state, DIMS(state)
     &     )

      implicit none

#include "REAL.H"
#include "HEATTRANSFER_F.H"
#include "ChemDriver_F.H"
#include "cdwrk.H"

c     arguments
      integer ncomps
      integer nspecies
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer DIMDEC(state)
      integer index_of_firstY
      integer index_of_lastY
      integer index_of_rho
      integer index_of_rhoH
      integer index_of_T
      integer maximum_iterations
      REAL_T maximum_error
      REAL_T dx(BL_SPACEDIM)
      REAL_T state(DIMV(state),ncomps)

c     local variables
      external FORT_TfromHYpt 
      integer FORT_TfromHYpt 
      integer i
      integer iterations
      integer j
      integer n
      integer offset
      logical first
      REAL_T ambient_pressure_cgs
      REAL_T ambient_pressure_si
      REAL_T contiguous_X
      REAL_T contiguous_Y
      REAL_T cpms
      REAL_T D
      REAL_T dpdt_factor
      REAL_T floored_T
      REAL_T FORT_P1ATMMKS
      REAL_T H
      REAL_T mmw
      REAL_T one_atmosphere_cgs
      REAL_T one_atmosphere_si
      REAL_T one_over_rho
      REAL_T one_over_wt
      REAL_T residual
      REAL_T rho_cgs
      REAL_T RU
      REAL_T RUC
      REAL_T smallest_time_step
      REAL_T space_step
      REAL_T T
      REAL_T time_step
      REAL_T Upsilon
      REAL_T wt
      allocatable contiguous_X
      allocatable contiguous_Y
      allocatable cpms
      allocatable one_over_wt
      allocatable residual
      allocatable Upsilon
      allocatable wt
      dimension contiguous_X (:)
      dimension contiguous_Y (:)
      dimension cpms (:)
      dimension one_over_wt (:)
      dimension residual (:)
      dimension Upsilon (:)
      dimension wt (:)

      allocate (contiguous_X (nspecies))
      allocate (contiguous_Y (nspecies))
      allocate (cpms (nspecies))
      allocate (one_over_wt (nspecies))
      allocate (residual (maximum_iterations))
      allocate (Upsilon (nspecies))
      allocate (wt (nspecies))

c///
c///  initialization
c///

c     offset of species in the state
      offset = index_of_firstY - 1

c     get ambient pressure in cgs units
c     one_atmosphere_cgs = 1013250.0 (dyne / cm**2)
c     one_atmosphere_si = 101325.0 (pascal)
      call CKRP (IWRK(ckbi), RWRK(ckbr), RU, RUC, one_atmosphere_cgs)
      call FORT_GETPAMB (ambient_pressure_si, dpdt_factor)
      one_atmosphere_si = FORT_P1ATMMKS ()
      ambient_pressure_cgs = one_atmosphere_cgs * ambient_pressure_si / one_atmosphere_si

c     reciprocal of molecular weights
      call CKWT (IWRK(ckbi),RWRK(ckbr),wt)
      do n=1,nspecies
         one_over_wt(n) = 1.0d0 / wt(n)
      end do

c     smallest space step
      space_step = min (dx(1), dx(2))

c///
c///  find the smallest time step for the given part of the domain
c///

      first = .true.
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

c           calculate T and Y from rho, rho H, and rho Y in the state
            one_over_rho = 1.0d0 / state(i, j, index_of_rho)
            H = one_over_rho * state(i, j, index_of_rhoH)
            do n=1,nspecies
               contiguous_Y(n) = one_over_rho * state(i, j, offset + n)
            end do
c           use the state value for T as a guess for the Newton iteration
            T = state(i, j, index_of_T)
            iterations = FORT_TfromHYpt 
     &           (T, 
     &           H, 
     &           contiguous_Y,
     &           maximum_error,
     &           maximum_iterations,
     &           residual)
            floored_T = MAX (T, TMIN_TRANS) 

c           evaluate the diffusion coefficients
            CALL CKCPMS 
     &           (floored_T, IWRK(ckbi), RWRK(ckbr), cpms)
            CALL CKYTX 
     &           (contiguous_Y, IWRK(ckbi), RWRK(ckbr), contiguous_X)
            CALL EGSPAR 
     &           (floored_T, contiguous_X, contiguous_Y, cpms, RWRK(egbr), 
     &           IWRK(egbi))
            CALL CKMMWY 
     &           (contiguous_Y, IWRK(ckbi), RWRK(ckbr), mmw)
            CALL EGSV1 
     &           (ambient_pressure_cgs, floored_T, contiguous_Y, mmw, 
     &           RWRK(egbr), Upsilon)
c           get a rho_cgs consistent with floored_T and the rest
            CALL CKRHOY 
     &           (ambient_pressure_cgs, floored_T, contiguous_Y, IWRK(ckbi),
     &           RWRK(ckbr), rho_cgs)

c           select the smallest time step
            do n=1,nspecies
c              per page 17 in the EGLIB manual, Upsilon has units (cm**2 / s).  
c              per equation 12 in the EGLIB manual, the mixture averaged diffusion 
c              coefficient is Upsilon times the dimensionless mmw / wt(n).  thus 
c              the resulting number with units (cm**2 / s) = 1.0d-4 (m**2 / s) must
c              be multiplied by 1.0d-4 to obtain si units.
               D = mmw * one_over_wt(n) * Upsilon(n) * 1.0d-4
               time_step = space_step ** 2 / (D * 2.0d0 * BL_SPACEDIM)
               if (first) then
                  smallest_time_step = time_step
                  first = .false.
               else
                  smallest_time_step = min (time_step, smallest_time_step)
               end if
            end do

         end do
      end do

c///
c///  release temporary arrays (not strictly necessary)
c///

      deallocate (contiguous_x)
      deallocate (contiguous_Y)
      deallocate (cpms)
      deallocate (one_over_wt)
      deallocate (residual)
      deallocate (Upsilon)
      deallocate (wt)

c///
c///  end of FORT_RK_STEP_SELECTION
c///

      return
      end
