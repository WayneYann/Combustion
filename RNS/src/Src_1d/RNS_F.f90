
subroutine rns_dudt_ad (lo, hi, &
     U, U_l1, U_h1, &
     dUdt, Ut_l1, Ut_h1, &
     flux, f_l1, f_h1, &
     dx)
  use meth_params_module, only : NVAR, gravity, URHO, UMX, UEDEN
  use hypterm_module, only : hypterm
  use difterm_module, only : difterm
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  U_l1,  U_h1
  integer, intent(in) :: Ut_l1, Ut_h1
  integer, intent(in) ::  f_l1,  f_h1
  double precision, intent(in)    ::    U( U_l1: U_h1,NVAR)
  double precision, intent(inout) :: dUdt(Ut_l1:Ut_h1,NVAR)
  double precision, intent(  out) :: flux( f_l1: f_h1,NVAR)
  double precision, intent(in) :: dx(1)

  integer :: Ulo(1), Uhi(1), i, n
  double precision :: dxinv(1)
  double precision, allocatable :: fdif(:,:)
  
  dxinv(1) = 1.d0/dx(1)
  
  Ulo(1) = U_l1
  Uhi(1) = U_h1

  allocate(fdif(lo(1):hi(1)+1,NVAR))
  
  if (f_l1.ne.lo(1) .or. f_h1.ne.hi(1)+1) then
     print *, 'flux has wrong size!'
     stop
  end if

  call hypterm(lo,hi,U,Ulo,Uhi,flux, dx)
  call difterm(lo,hi,U,Ulo,Uhi,fdif, dxinv)
  
  do n=1, NVAR
     flux(lo(1),n) = flux(lo(1),n) + fdif(lo(1),n)
     do i=lo(1),hi(1)
        flux(i+1,n) = flux(i+1,n) + fdif(i+1,n)
        dUdt(i,n) = dxinv(1) * (flux(i,n) - flux(i+1,n))
     end do
  end do
  
  deallocate(fdif)

  if (gravity .ne. 0.d0) then
     do i=lo(1),hi(1)
        dUdt(i,UMX  ) = dUdt(i,UMX  ) + U(i,URHO)*gravity
        dUdt(i,UEDEN) = dUdt(i,UEDEN) + U(i,UMX )*gravity
     end do
  end if

end subroutine rns_dudt_ad

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_advchem(lo,hi,U,U_l1,U_h1,dt)
  use meth_params_module, only : NVAR
  use chemterm_module, only : chemterm
  implicit none
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  U_l1,  U_h1
  double precision, intent(inout) :: U(U_l1: U_h1,NVAR)
  double precision, intent(in) :: dt

  integer :: Ulo(1), Uhi(1)

  Ulo(1) = U_l1
  Uhi(1) = U_h1
  call chemterm(lo, hi, U, Ulo, Uhi, dt)
end subroutine rns_advchem

subroutine rns_advchem2(lo,hi,U,U_l1,U_h1,Up,Up_l1,Up_h1,dt)
  use meth_params_module, only : NVAR
  use chemterm_module, only : chemterm
  implicit none
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  U_l1, U_h1, Up_l1, Up_h1
  double precision, intent(inout) :: U(U_l1: U_h1,NVAR)
  double precision, intent(in) :: Up(Up_l1: Up_h1,NVAR)
  double precision, intent(in) :: dt

  integer :: Ulo(1), Uhi(1)

  Ulo(1) = U_l1
  Uhi(1) = U_h1
  call chemterm(lo, hi, U, Ulo, Uhi, dt, Up)
end subroutine rns_advchem2

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_dUdt_chem(lo,hi,U,U_l1,U_h1,Ut,Ut_l1,Ut_h1)
  use meth_params_module, only : NVAR
  use chemterm_module, only : dUdt_chem
  implicit none
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  U_l1,  U_h1, Ut_l1,  Ut_h1
  double precision, intent(in ) ::  U( U_l1: U_h1,NVAR)
  double precision, intent(out) :: Ut(Ut_l1:Ut_h1,NVAR)

  integer :: Ulo(1), Uhi(1), Utlo(1), Uthi(1)

  Ulo(1) = U_l1
  Uhi(1) = U_h1
  Utlo(1) = Ut_l1
  Uthi(1) = Ut_h1
  call dUdt_chem(lo, hi, U, Ulo, Uhi, Ut, Utlo, Uthi)
end subroutine rns_dUdt_chem

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

  integer :: i, pt_index(1), ierr
  double precision :: rhoInv, e, v, Y(NSPEC)

  do i=lo(1),hi(1)
     pt_index(1) = i

     rhoInv = 1.0d0/U(i,URHO)

     v  = U(i,UMX)*rhoInv     
     e  = U(i,UEDEN)*rhoInv - 0.5d0*v*v

     Y = U(i,UFS:UFS+NSPEC-1)*rhoInv

     call eos_get_T(U(i,UTEMP), e, Y, pt_index, ierr)

     if (ierr .ne. 0) then
        print *, 'rns_compute_temp failed at ', i,U(i,:)
        call bl_error("rns_compute_temp failed")
     end if
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
  double precision :: dom_spec,x,rhoInv, sumrY, fac

  double precision, parameter :: eps = -1.d-16

  do i = lo(1),hi(1)

     any_negative = .false.

     rhoInv = 1.d0/U(i,URHO)

     sumrY = 0.d0

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

        sumrY = sumrY + U(i,n)
     end do

     fac = U(i,URHO)/sumrY
     do n = UFS, UFS+nspec-1
        U(i,n) = U(i,n)*fac
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
              
              x = U(i,n)*rhoInv
              
              ! ! Here we only print the bigger negative values
              ! if (x .lt. -1.d-2) then
              !    print *,'Correcting negative species   ',n-UFS+1
              !    print *,'   at cell (i)                ',i
              !    print *,'Negative (rho*Y) is           ',U(i,n)
              !    print *,'Negative      Y  is           ',x
              !    print *,'Filling from dominant species ',int_dom_spec-UFS+1
              !    print *,'  which had Y =               ',&
              !         U(i,int_dom_spec) / U(i,URHO)
              ! end if

              ! Take enough from the dominant species to fill the negative one.
              U(i,int_dom_spec) = U(i,int_dom_spec) + U(i,n)
   
              ! Test that we didn't make the dominant species negative
              if (U(i,int_dom_spec) .lt. 0.d0) then 
                 print *,' Just made dominant species negative ',int_dom_spec+UFS-1,' at ',i
                 print *,'We were fixing species ',n-UFS+1,' which had value ',x
                 print *,'Dominant species became ',U(i,int_dom_spec)*rhoInv
                 call bl_error("Error:: CNSReact_2d.f90 :: ca_enforce_nonnegative_species")
              end if

              ! Now set the negative species to zero
              U(i,n) = 0.d0

           end if

        end do

     end if
     
  end do

end subroutine rns_enforce_consistent_Y


subroutine rns_sum_cons ( &
     U  ,U_l1,U_h1, &
     msk,m_l1,m_h1, &
     vol,v_l1,v_h1, &
     s)
  use meth_params_module, only : NVAR
  implicit none
  
  integer, intent(in) :: U_l1,U_h1
  integer, intent(in) :: m_l1,m_h1
  integer, intent(in) :: v_l1,v_h1
  double precision, intent(in) :: U  (U_l1:U_h1,NVAR)
  double precision, intent(in) :: msk(m_l1:m_h1)
  double precision, intent(in) :: vol(v_l1:v_h1)
  double precision, intent(inout) :: s(3)

  integer :: i, n

  do n=1,3
     do i=m_l1,m_h1
        s(n) = s(n) + msk(i)*vol(i)*U(i,n)
     end do
  end do

end subroutine rns_sum_cons


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
