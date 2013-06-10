! :: ----------------------------------------------------------
! :: Volume-weight average the fine grid data onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  ngc        => number of ghost cells in coarse array
! ::  nvar	 => number of components in arrays
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  ngf        => number of ghost cells in fine array
! ::  rfine      => (ignore) used in 2-D RZ calc
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
      subroutine rns_avgdown (crse,c_l1,c_h1,nvar, &
           cv,cv_l1,cv_h1, &
           fine,f_l1,f_h1, &
           fv,fv_l1,fv_h1,lo,hi,lrat)

      implicit none
      integer c_l1,c_h1
      integer cv_l1,cv_h1
      integer f_l1,f_h1
      integer fv_l1,fv_h1
      integer lo(1), hi(1)
      integer nvar, lrat(1)
      double precision crse(c_l1:c_h1,nvar)
      double precision cv(cv_l1:cv_h1)
      double precision fine(f_l1:f_h1,nvar)
      double precision fv(fv_l1:fv_h1)

      integer i, n, ic, ioff
      integer lratx
 
      lratx = lrat(1)
 
      do n = 1, nvar
 
!        Set coarse grid to zero on overlap
         do ic = lo(1), hi(1)
            crse(ic,n) = 0.d0
         enddo
 
 
!        Sum fine data
         do ioff = 0, lratx-1
            do ic = lo(1), hi(1)
               i = ic*lratx + ioff
               crse(ic,n) = crse(ic,n) + fv(i) * fine(i,n)
            enddo
         enddo
             
!        Divide out by volume weight
         do ic = lo(1), hi(1)
            crse(ic,n) = crse(ic,n) / cv(ic)
         enddo
            
      enddo

      end subroutine rns_avgdown


! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine rns_enforce_nonnegative_species(uout,uout_l1,uout_h1,lo,hi)

      end subroutine rns_enforce_nonnegative_species

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine rns_estdt(u,u_l1,u_h1,lo,hi,dx,dt)
        use eos_module
        use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC
        implicit none

        integer u_l1,u_h1
        integer lo(1), hi(1)
        double precision u(u_l1:u_h1,NVAR)
        double precision dx(1), dt

        integer :: i
        double precision :: rhoInv, ux, T, e, c, xn(NSPEC)

        do i = lo(1), hi(1)
           rhoInv = 1.d0/u(i,URHO)

         ux = u(i,UMX)*rhoInv
         T  = u(i,UTEMP)

         e = u(i,UEDEN)*rhoInv - 0.5d0*ux*ux

         if (NSPEC > 0) then
            xn = u(i,UFS:UFS+NSPEC-1)*rhoInv
         end if

         call eos_get_soundspeed(c,u(i,URHO),e,T,xn)

        end do

      end subroutine rns_estdt
