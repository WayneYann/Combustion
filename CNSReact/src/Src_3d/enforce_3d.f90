! ::
! :: ----------------------------------------------------------
! ::

subroutine cns_enforce_nonnegative_species( u, &
     u_l1,u_l2,u_l3, &
     u_h1,u_h2,u_h3,lo,hi)

  use chemistry_module, only : nspecies
  use meth_params_module, only : NVAR, URHO, UFS

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: u_l1, u_l2, u_l3, u_h1, u_h2, u_h3
  double precision, intent(inout) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)

  ! Local variables
  integer          :: i,j,k,n
  integer          :: int_dom_spec
  logical          :: any_negative
  double precision :: dom_spec,x

  double precision, parameter :: eps = -1.0d-16

  !$OMP PARALLEL DO PRIVATE(i,j,k,n,dom_spec,int_dom_spec,any_negative,x)
  do k = lo(3),hi(3)
  do j = lo(2),hi(2)
  do i = lo(1),hi(1)

     any_negative = .false.
     !
     ! First deal with tiny undershoots by just setting them to zero.
     !
     do n = UFS, UFS+nspecies-1
        if (u(i,j,k,n) .lt. 0.d0) then
           x = u(i,j,k,n)/u(i,j,k,URHO)
           if (x .gt. eps) then
              u(i,j,k,n) = 0.d0
           else
              any_negative = .true.
           end if
        end if
     end do
     !
     ! We know there are one or more undershoots needing correction.
     !
     if (any_negative) then
        !
        ! Find the dominant species.
        !
        int_dom_spec = UFS
        dom_spec     = u(i,j,k,int_dom_spec)

        do n = UFS,UFS+nspecies-1
           if (u(i,j,k,n) .gt. dom_spec) then
              dom_spec     = u(i,j,k,n)
              int_dom_spec = n
           end if
        end do
        !
        ! Now take care of undershoots greater in magnitude than 1e-16.
        !
        do n = UFS, UFS+nspecies-1

           if (u(i,j,k,n) .lt. 0.d0) then

              x = u(i,j,k,n)/u(i,j,k,URHO)
              !
              ! Here we only print the bigger negative values.
              !
              if (x .lt. -1.d-2) then
                 print *,'Correcting negative species   ',n
                 print *,'   at cell (i,j,k)            ',i,j,k
                 print *,'Negative (rho*X) is           ',u(i,j,k,n)
                 print *,'Negative      X  is           ',x
                 print *,'Filling from dominant species ',int_dom_spec
                 print *,'  which had X =               ',&
                      u(i,j,k,int_dom_spec) / u(i,j,k,URHO)
              end if
              !
              ! Take enough from the dominant species to fill the negative one.
              !
              u(i,j,k,int_dom_spec) = u(i,j,k,int_dom_spec) + u(i,j,k,n)
              !
              ! Test that we didn't make the dominant species negative.
              !
              if (u(i,j,k,int_dom_spec) .lt. 0.d0) then 
                 print *,' Just made dominant species negative ',int_dom_spec,' at ',i,j,k 
                 print *,'We were fixing species ',n,' which had value ',x
                 print *,'Dominant species became ',u(i,j,k,int_dom_spec) / u(i,j,k,URHO)
                 call bl_error("Error:: CNSReact_3d.f90 :: cns_enforce_nonnegative_species")
              end if
              !
              ! Now set the negative species to zero.
              !
              u(i,j,k,n) = 0.d0

           end if

        enddo
     end if
  enddo
  enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine cns_enforce_nonnegative_species

!
! reset total energy
!
subroutine cns_enforce_consistent_e(lo,hi,S, &
     S_l1,S_l2,S_l3,S_h1,S_h2,S_h3)

  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: S_l1,S_l2,S_l3,S_h1,S_h2,S_h3
  double precision, intent(inout) :: S(S_l1:S_h1,S_l2:S_h2,S_l3:S_h3,NVAR)

  ! Local variables
  integer          :: i,j,k

  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           S(i,j,k,UEDEN) = S(i,j,k,UEINT) + &
                0.5d0*(S(i,j,k,UMX)**2 &
                &     +S(i,j,k,UMY)**2 &
                &     +S(i,j,k,UMZ)**2) / S(i,j,k,URHO)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine cns_enforce_consistent_e

! :::
! ::: ------------------------------------------------------------------
! :::
subroutine cns_reset_internal_energy(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi,verbose)

  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT

  implicit none

  integer, intent(in) :: lo(3), hi(3), verbose
  integer, intent(in) :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
  double precision, intent(inout) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)

  integer          :: i,j,k

  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           u(i,j,k,UEINT) = u(i,j,k,UEDEN) - &
                0.5d0*(u(i,j,k,UMX)**2 &
                &     +u(i,j,k,UMY)**2 &
                &     +u(i,j,k,UMZ)**2) / u(i,j,k,URHO)
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine cns_reset_internal_energy
