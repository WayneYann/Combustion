
subroutine rns_dudt (lo, hi, &
     U, U_l1, U_h1, &
     dUdt, Ut_l1, Ut_h1, &
     dx, dt)
  use meth_params_module, only : NVAR
  use weno_module, only : reconstruct
  use riemann_module, only : riemann
  use chemterm_module, only : chemterm
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  U_l1,  U_h1
  integer, intent(in) :: Ut_l1, Ut_h1
  double precision, intent(in)    ::    U( U_l1: U_h1,NVAR)
  double precision, intent(inout) :: dUdt(Ut_l1:Ut_h1,NVAR)
  double precision, intent(in) :: dx(1), dt

  integer :: Ulo(1), Uhi(1), Utlo(1), Uthi(1), i, n
  double precision :: dxinv(1)
  double precision, allocatable :: UL(:,:), UR(:,:), fx(:,:)

  dxinv(1) = 1.d0/dx(1)

  Ulo(1) = U_l1
  Uhi(1) = U_h1

  Utlo(1) = Ut_l1
  Uthi(1) = Ut_h1

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

  call chemterm(lo, hi, U, Ulo, Uhi, dUdt, Utlo, Uthi, dt)

  deallocate(UL,UR,fx)

end subroutine rns_dudt

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_compute_temp(lo,hi,U,U_l1,U_h1)
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC
  use eos_module, only : eos_get_T
  implicit none
  
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  U_l1,  U_h1
  double precision, intent(inout) :: U( U_l1: U_h1,NVAR)

  integer :: i
  double precision :: rhoInv, e, v, Y(NSPEC)

  do i=lo(1),hi(1)
     rhoInv = 1.0d0/U(i,URHO)

     v  = U(i,UMX)*rhoInv     
     e  = U(i,UEDEN)*rhoInv - 0.5d0*v*v

     Y = U(i,UFS:UFS+NSPEC-1)*rhoInv

     call eos_get_T(U(i,UTEMP), e, Y)
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

  ! Local variables
  integer          :: i,n
  integer          :: int_dom_spec
  logical          :: any_negative
  double precision :: dom_spec,x,rhoInv

  double precision, parameter :: eps = -1.d-16

  do i = lo(1),hi(1)

     any_negative = .false.

     rhoInv = 1.d0/U(i,URHO)

     ! First deal with tiny undershoots by just setting them to zero
     do n = UFS, UFS+nspec-1
        if (U(i,n) .lt. 0.d0) then
           x = U(i,n) * rhoInv
           if (x .gt. eps) then
              U(i,n) = 0.d0
           else
              any_negative = .true.
           end if
        end if
     end do

     ! We know there are one or more undershoots needing correction 
     if (any_negative) then

        ! Find the dominant species
        dom_spec = 0.d0
        int_dom_spec = 0
        do n = UFS,UFS+nspec-1
           if (U(i,n) .gt. dom_spec) then
              dom_spec = U(i,n)
              int_dom_spec = n
           end if
        end do

        ! Now take care of undershoots greater in magnitude than 1e-16.
        do n = UFS, UFS+nspec-1
           
           if (U(i,n) .lt. 0.d0) then
              
              x = U(i,n)/U(i,URHO)
              
              ! Here we only print the bigger negative values
              if (x .lt. -1.d-2) then
                 print *,'Correcting negative species   ',n
                 print *,'   at cell (i)                ',i
                 print *,'Negative (rho*X) is           ',U(i,n)
                 print *,'Negative      X  is           ',x
                 print *,'Filling from dominant species ',int_dom_spec
                 print *,'  which had X =               ',&
                      U(i,int_dom_spec) / U(i,URHO)
              end if

              ! Take enough from the dominant species to fill the negative one.
              U(i,int_dom_spec) = U(i,int_dom_spec) + U(i,n)
   
              ! Test that we didn't make the dominant species negative
              if (U(i,int_dom_spec) .lt. 0.d0) then 
                 print *,' Just made dominant species negative ',int_dom_spec,' at ',i
                 print *,'We were fixing species ',n,' which had value ',x
                 print *,'Dominant species became ',U(i,int_dom_spec) / U(i,URHO)
                 call bl_error("Error:: CNSReact_2d.f90 :: ca_enforce_nonnegative_species")
              end if

              ! Now set the negative species to zero
              U(i,n) = 0.d0

           end if

        end do

     end if
     
  end do

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
        double precision :: rhoInv, ux, T, e, c, Y(NSPEC)

        do i = lo(1), hi(1)
           rhoInv = 1.d0/u(i,URHO)

           ux = u(i,UMX)*rhoInv
           T  = u(i,UTEMP)
           
           e = u(i,UEDEN)*rhoInv - 0.5d0*ux*ux
           
           if (NSPEC > 0) then
              Y = u(i,UFS:UFS+NSPEC-1)*rhoInv
           end if
           
           call eos_get_c(c,u(i,URHO),T,Y)
           
           dt = min(dt, dx(1)/(abs(ux)+c+1.d-50))
        end do

      end subroutine rns_estdt
