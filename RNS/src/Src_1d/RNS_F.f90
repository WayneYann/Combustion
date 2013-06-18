
subroutine rns_dudt (lo, hi, &
     U, U_l1, U_h1, &
     dUdt, Ut_l1, Ut_h1, &
     dx)
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC
  use weno_module, only : reconstruct
  use riemann_module, only : riemann
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  U_l1,  U_h1
  integer, intent(in) :: Ut_l1, Ut_h1
  double precision, intent(in)    ::    U( U_l1: U_h1,NVAR)
  double precision, intent(inout) :: dUdt(Ut_l1:Ut_h1,NVAR)
  double precision, intent(in) :: dx(1)

  integer :: Ulo(1), Uhi(1), i, n
  double precision :: dxinv(1)
  double precision, allocatable :: UL(:,:), UR(:,:), fx(:,:)

  dxinv(1) = 1.d0/dx(1)

  Ulo(1) = U_l1
  Uhi(1) = U_h1

  allocate(UL(lo(1):hi(1)+1,NVAR))
  allocate(UR(lo(1):hi(1)+1,NVAR))
  allocate(fx(lo(1):hi(1)+1,NVAR))

  call reconstruct(lo, hi, U, Ulo, Uhi, UL, UR)

  call riemann(lo, hi, UL, UR, fx)

  do n=1, NVAR
     do i=lo(1),hi(1)
        dUdt(i,n) = dxinv(1) * (fx(i,n) - fx(i+1,n))
     end do
  end do

  deallocate(UL,UR,fx)

end subroutine rns_dudt

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_compute_temp(lo,hi,U,U_l1,U_h1)
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC
  use eos_module, only : eos_given_ReY
  implicit none
  
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  U_l1,  U_h1
  double precision, intent(inout) :: U( U_l1: U_h1,NVAR)

  integer :: i
  double precision :: gamc, p, c, rho, rhoInv, e, v, Y(NSPEC)

  do i=lo(1),hi(1)

     rho  = U(i,URHO)
     rhoInv = 1.0d0/rho

     v  = U(i,UMX)*rhoInv
     
     e  = U(i,UEDEN)*rhoInv - 0.5d0*v*v

     if (NSPEC > 0) then
        Y = U(i,UFS:UFS+NSPEC-1)*rhoInv
     end if

     call eos_given_ReY(gamc,p,c, U(i,UTEMP), rho, e, Y)

  end do
end subroutine rns_compute_temp

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_enforce_consistent_Y(lo,hi,U,U_l1,U_h1)
  use meth_params_module, only : NVAR, URHO, UFS, NSPEC
  implicit none
  
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  U_l1,  U_h1
  double precision, intent(inout) :: U( U_l1: U_h1,NVAR)
  print *, 'rns_enforce_consistent_Y not implemented'
  stop
end subroutine rns_enforce_consistent_Y



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

      subroutine rns_estdt(u,u_l1,u_h1,lo,hi,dx,dt)
        use eos_module, only : eos_get_c
        use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC
        implicit none

        integer u_l1,u_h1
        integer lo(1), hi(1)
        double precision u(u_l1:u_h1,NVAR)
        double precision dx(1), dt

        integer :: i
        double precision :: rhoInv, ux, T, e, p, c, g, Y(NSPEC)

        do i = lo(1), hi(1)
           rhoInv = 1.d0/u(i,URHO)

           ux = u(i,UMX)*rhoInv
           T  = u(i,UTEMP)
           
           e = u(i,UEDEN)*rhoInv - 0.5d0*ux*ux
           
           if (NSPEC > 0) then
              Y = u(i,UFS:UFS+NSPEC-1)*rhoInv
           end if
           
           call eos_get_c(c,u(i,URHO),T,Y)
           
           dt = min(dt, dx(1)/(c+1.d-50))
        end do

      end subroutine rns_estdt
