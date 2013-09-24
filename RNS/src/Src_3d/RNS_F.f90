
subroutine rns_dudt (lo, hi, &
     U   , U_l1, U_l2, U_l3, U_h1, U_h2, U_h3, &
     dUdt,Ut_l1,Ut_l2,Ut_l3,Ut_h1,Ut_h2,Ut_h3, &
     xflx,xf_l1,xf_l2,xf_l3,xf_h1,xf_h2,xf_h3, &
     yflx,yf_l1,yf_l2,yf_l3,yf_h1,yf_h2,yf_h3, &
     zflx,zf_l1,zf_l2,zf_l3,zf_h1,zf_h2,zf_h3, &
     dx)
  use meth_params_module, only : NVAR, gravity, URHO, UMZ, UEDEN, &
       xblksize, yblksize, zblksize, nthreads
  use hypterm_module, only : hypterm
  use difterm_module, only : difterm
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) ::  U_l1,  U_h1,  U_l2,  U_h2,  U_l3,  U_h3
  integer, intent(in) :: Ut_l1, Ut_h1, Ut_l2, Ut_h2, Ut_l3, Ut_h3
  integer, intent(in) :: xf_l1, xf_h1, xf_l2, xf_h2, xf_l3, xf_h3
  integer, intent(in) :: yf_l1, yf_h1, yf_l2, yf_h2, yf_l3, yf_h3
  integer, intent(in) :: zf_l1, zf_h1, zf_l2, zf_h2, zf_l3, zf_h3
  double precision,intent(in) ::   U( U_l1: U_h1, U_l2: U_h2, U_l3: U_h3,NVAR)
  double precision,intent(out)::dUdt(Ut_l1:Ut_h1,Ut_l2:Ut_h2,Ut_l3:Ut_h3,NVAR)
  double precision,intent(out)::xflx(xf_l1:xf_h1,xf_l2:xf_h2,xf_l3:xf_h3,NVAR)
  double precision,intent(out)::yflx(yf_l1:yf_h1,yf_l2:yf_h2,yf_l3:yf_h3,NVAR)
  double precision,intent(out)::zflx(zf_l1:zf_h1,zf_l2:zf_h2,zf_l3:zf_h3,NVAR)
  double precision,intent(in) :: dx(3)

  integer :: Ulo(3), Uhi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), tlo(3), thi(3)
  integer :: i, j, k, n, blocksize(3), ib, jb, kb, nb(3)
  double precision :: dxinv(3)
  double precision, allocatable :: bxflx(:,:,:,:), byflx(:,:,:,:), bzflx(:,:,:,:)

  xflx = 0.d0
  yflx = 0.d0
  zflx = 0.d0

  dUdt = 0.d0

end subroutine rns_dudt

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_advchem(lo,hi,U,U_l1,U_l2,U_l3,U_h1,U_h2,U_h3,dt)
  use meth_params_module, only : NVAR
  use chemterm_module, only : chemterm
  implicit none
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) ::  U_l1, U_l2, U_l3, U_h1, U_h2, U_h3
  double precision, intent(inout) :: U(U_l1:U_h1,U_l2:U_h2,U_l3:U_h3,NVAR)
  double precision, intent(in) :: dt

  integer :: Ulo(3), Uhi(3)

  Ulo(1) = U_l1
  Ulo(2) = U_l2
  Ulo(3) = U_l3
  Uhi(1) = U_h1
  Uhi(2) = U_h2
  Uhi(3) = U_h3
  call chemterm(lo, hi, U, Ulo, Uhi, dt)
end subroutine rns_advchem

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_compute_temp(lo,hi,U,U_l1,U_l2,U_l3,U_h1,U_h2,U_h3)
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, NSPEC
  use eos_module, only : eos_get_T
  implicit none
  
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: U_l1, U_l2, U_l3, U_h1, U_h2, U_h3
  double precision, intent(inout) :: U(U_l1:U_h1,U_l2:U_h2,U_l3:U_h3,NVAR)

  integer :: i, j, k
  double precision :: rhoInv, e, vx, vy, vz, Y(NSPEC)

  !$omp parallel do private(i,j,k,rhoInv,e,vx,vy,vz,Y) collapse(2)
  do k=lo(3),hi(3)
  do j=lo(2),hi(2)
  do i=lo(1),hi(1)
     rhoInv = 1.0d0/U(i,j,k,URHO)

     vx = U(i,j,k,UMX)*rhoInv     
     vy = U(i,j,k,UMY)*rhoInv     
     vz = U(i,j,k,UMZ)*rhoInv     
     e  = U(i,j,k,UEDEN)*rhoInv - 0.5d0*(vx**2+vy**2+vz**2)

     Y = U(i,j,k,UFS:UFS+NSPEC-1)*rhoInv

     call eos_get_T(U(i,j,k,UTEMP), e, Y)
  end do
  end do
  end do
  !$omp end parallel do
end subroutine rns_compute_temp

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine rns_enforce_consistent_Y(lo,hi,U,U_l1,U_l2,U_l3,U_h1,U_h2,U_h3)
  use meth_params_module, only : NVAR, URHO, UFS, NSPEC
  implicit none
  
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: U_l1, U_l2, U_l3, U_h1, U_h2, U_h3
  double precision, intent(inout) :: U(U_l1:U_h1,U_l2:U_h2,U_l3:U_h3,NVAR)

  ! Local variables
  integer          :: i,j,k,n
  integer          :: int_dom_spec
  logical          :: any_negative
  double precision :: dom_spec,x,rhoInv, sumrY, fac

  double precision, parameter :: eps = -1.d-16

  !$omp parallel do private(i,j,k,n,int_dom_spec,any_negative,dom_spec) &
  !$omp private(x,rhoInv,sumrY,fac) collapse(2)
  do k = lo(3),hi(3)
  do j = lo(2),hi(2)
  do i = lo(1),hi(1)

     any_negative = .false.

     rhoInv = 1.d0/U(i,j,k,URHO)

     sumrY = 0.d0

     ! First deal with tiny undershoots by just setting them to zero
     do n = UFS, UFS+nspec-1
        if (U(i,j,k,n) .lt. 0.d0) then
           x = U(i,j,k,n) * rhoInv
           if (x .gt. eps) then
              U(i,j,k,n) = 0.d0
           else
              any_negative = .true.
           end if
        end if

        sumrY = sumrY + U(i,j,k,n)
     end do

     fac = U(i,j,k,URHO)/sumrY
     do n = UFS, UFS+nspec-1
        U(i,j,k,n) = U(i,j,k,n)*fac
     end do

     ! We know there are one or more undershoots needing correction 
     if (any_negative) then

        ! Find the dominant species
        dom_spec = 0.d0
        int_dom_spec = 0
        do n = UFS,UFS+nspec-1
           if (U(i,j,k,n) .gt. dom_spec) then
              dom_spec = U(i,j,k,n)
              int_dom_spec = n
           end if
        end do

        ! Now take care of undershoots greater in magnitude than 1e-16.
        do n = UFS, UFS+nspec-1
           
           if (U(i,j,k,n) .lt. 0.d0) then
              
              x = U(i,j,k,n)*rhoInv
              
              ! ! Here we only print the bigger negative values
              ! if (x .lt. -1.d-2) then
              !    print *,'Correcting negative species   ',n-UFS+1
              !    print *,'   at cell (i,j,k)            ',i,j
              !    print *,'Negative (rho*Y) is           ',U(i,j,k,n)
              !    print *,'Negative      Y  is           ',x
              !    print *,'Filling from dominant species ',int_dom_spec-UFS+1
              !    print *,'  which had Y =               ',&
              !         U(i,j,k,int_dom_spec) / U(i,j,k,URHO)
              ! end if

              ! Take enough from the dominant species to fill the negative one.
              U(i,j,k,int_dom_spec) = U(i,j,k,int_dom_spec) + U(i,j,k,n)
   
              ! ! Test that we didn't make the dominant species negative
              ! if (U(i,j,k,int_dom_spec) .lt. 0.d0) then 
              !    print *,' Just made dominant species negative ',int_dom_spec-UFS+1,' at ',i
              !    print *,'We were fixing species ',n-UFS+1,' which had value ',x
              !    print *,'Dominant species became ',U(i,j,k,int_dom_spec) / U(i,j,k,URHO)
              !    call bl_error("Error:: CNSReact_2d.f90 :: ca_enforce_nonnegative_species")
              ! end if

              ! Now set the negative species to zero
              U(i,j,k,n) = 0.d0

           end if

        end do

     end if
     
  end do
  end do
  end do
  !$omp end parallel do

end subroutine rns_enforce_consistent_Y


! :: ----------------------------------------------------------
! :: Volume-weight average the fine grid data onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  nvar	 => number of components in arrays
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
      subroutine rns_avgdown(crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,nvar, &
                             cv,cv_l1,cv_l2,cv_l3,cv_h1,cv_h2,cv_h3, &
                             fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                             fv,fv_l1,fv_l2,fv_l3,fv_h1,fv_h2,fv_h3,lo,hi,lrat)

      implicit none
      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      integer cv_l1,cv_l2,cv_l3,cv_h1,cv_h2,cv_h3
      integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer fv_l1,fv_l2,fv_l3,fv_h1,fv_h2,fv_h3
      integer lo(3), hi(3)
      integer nvar, lrat(3)
      double precision crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,nvar)
      double precision cv(cv_l1:cv_h1,cv_l2:cv_h2,cv_l3:cv_h3)
      double precision fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3,nvar)
      double precision fv(fv_l1:fv_h1,fv_l2:fv_h2,fv_l3:fv_h3)

      integer i, j, k, n, ic, jc, kc, ioff, joff, koff
      integer lratx, lraty, lratz
      double precision   volfrac

      lratx   = lrat(1)
      lraty   = lrat(2)
      lratz   = lrat(3)
      volfrac = 1.d0/float(lrat(1)*lrat(2)*lrat(3))

      !$omp parallel do private(i,j,k,n,ic,jc,kc,ioff,joff,koff)
      do n = 1, nvar
         !
         ! Set coarse grid to zero on overlap.
         !
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = 0.d0
               enddo
            enddo
         enddo
         !
         ! Sum fine data.
         !
         do koff = 0, lratz-1
            do kc = lo(3),hi(3)
               k = kc*lratz + koff
               do joff = 0, lraty-1
                  do jc = lo(2), hi(2)
                     j = jc*lraty + joff
                     do ioff = 0, lratx-1
                        do ic = lo(1), hi(1)
                           i = ic*lratx + ioff
                           crse(ic,jc,kc,n) = crse(ic,jc,kc,n) + fine(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !
         ! Divide out by volume weight.
         !
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = volfrac*crse(ic,jc,kc,n)
               enddo
            enddo
         enddo

      enddo
      !$omp end parallel do

      end subroutine rns_avgdown

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine rns_estdt(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi,dx,dt)
        use eos_module, only : eos_get_c
        use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, NSPEC
        implicit none

        integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
        integer lo(3), hi(3)
        double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
        double precision dx(3), dt

        integer :: i, j, k
        double precision :: rhoInv, vx, vy, vz, T, e, c, Y(NSPEC)

        !$omp parallel do private(i,j,k,rhoInv,vx,vy,vz,T,e,c,Y) reduction(min:dt) &
        !$omp collapse(2)
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rhoInv = 1.d0/u(i,j,k,URHO)

           vx = u(i,j,k,UMX)*rhoInv
           vy = u(i,j,k,UMY)*rhoInv
           vz = u(i,j,k,UMZ)*rhoInv
           T  = u(i,j,k,UTEMP)
           
           e = u(i,j,k,UEDEN)*rhoInv - 0.5d0*(vx**2+vy**2+vz**2)
           
           Y = u(i,j,k,UFS:UFS+NSPEC-1)*rhoInv
           
           call eos_get_c(c,u(i,j,k,URHO),T,Y)

           dt = min(dt, dx(1)/(abs(vx)+c+1.d-50), &
                &       dx(2)/(abs(vy)+c+1.d-50), &
                &       dx(3)/(abs(vz)+c+1.d-50) )
        end do
        end do
        end do
        !$omp end parallel do

      end subroutine rns_estdt
